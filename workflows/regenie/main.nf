/** https://github.com/FINNGEN/regenie-pipelines */

include {SNPEFF_DOWNLOAD} from '../../modules/snpeff/download'



workflow {
	bgen = Channel.of(file(params.bgen),file(params.bgen.substring(0,params.bgen.lastIndexOf("."))+".sample")).collect()
	prune_ch = LD_PRUNING(bgen)
	vcf_ch = MAKE_VCF(bgen)
	snpeff_db = SNPEFF_DOWNLOAD(params.snpeff_db)
	STEP1(bgen,file(params.covarFile), file(params.phenoFile), prune_ch.output )

	ANNOT_VCF(vcf_ch.output, snpeff_db.output)
	}

/** STEP 1 of regenie cannot use more than 'X' alleles,  LD PRUNING is required */

process LD_PRUNING {
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(bgen_files)
output:
	path("keep.txt"),emit:output
script:
	def bgen = bgen_files.find{it.name.endsWith(".bgen")}
	def sample = bgen_files.find{it.name.endsWith(".sample")}
"""
mkdir -p TMP

set -x

# get rid of the sex; NA causes plink to fail ?
cut -d ' ' -f1,2,3 "${sample}" > TMP/jeter.sample

plink2 --bgen ${bgen} ref-first \\
	--sample TMP/jeter.sample \\
	--rm-dup force-first \\
    --allow-extra-chr \\
     --no-sex \\
    --indep-pairwise 50 1 0.2 \\
    --out TMP/plink


find TMP 1>&2

plink2 --bgen ${bgen} ref-first \\
	--sample TMP/jeter.sample \\
	 --allow-extra-chr  \\
	--extract TMP/plink.prune.in \\
	--make-bed \\
	--out TMP/pruned_data

plink --bfile TMP/pruned_data \\
	--allow-extra-chr \\
	--r2 --ld-window 50 \\
	--ld-window-kb 5000 \\
	--ld-window-r2 0.2 \\
	--out TMP/ld


awk '{print \$3;}'  TMP/ld.ld  | uniq  | sort -T . | uniq > TMP/snpInLD.txt


plink --bfile TMP/pruned_data \\
	--exclude TMP/snpInLD.txt \\
	--make-bed \\
	--out TMP/indepSNP_data

#
# it is not recommened to use more than 1000000 variants in step 1 of regenie
#

awk '{print \$2}' TMP/indepSNP_data.bim |\\
	LC_ALL=C sort -T TMP |\\
	head -n 1000000 > TMP/keep.txt

mv TMP/keep.txt ./

"""
}

process MAKE_VCF {
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(bgen_files)
output:
	path("sites.vcf.gz"),emit:output
script:
	        def bgen = bgen_files.find{it.name.endsWith(".bgen")}
        def sample = bgen_files.find{it.name.endsWith(".sample")}

"""
mkdir -p TMP
set -x

cut -d ' ' -f1,2,3 "${sample}" > TMP/jeter.sample

# 
# TODO : go faster and output only sites, no genotype ?
#
plink2 --bgen ${bgen} ref-first \\
        --sample TMP/jeter.sample \\
	--recode vcf bgz \\
        --out TMP/export

find TMP 1>&2
file TMP/export*

find TMP 1>&2
bcftools view -G -O u TMP/*.vcf.gz |\\
	bcftools annotate -x 'INFO,ID,QUAL,FILTER' -O z -o TMP/sites.vcf.gz

mv TMP/sites.vcf.gz ./

"""
}

process STEP1 {
conda "${moduleDir}/../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
	path(bgen_files)
        path(covariates)
	path(phenoFile)
	path(keep_rs)
output:
	path("${params.prefix}step1.log"),emit:log
script:
	def bgen = bgen_files.find{it.name.endsWith(".bgen")}
	def sample = bgen_files.find{it.name.endsWith(".sample")}
	def args = "--bsize 1000 --bt --phenoCol Y1 "
"""

mkdir -p TMP/OUT
set -x

cut -d ' ' -f1,2,3 "${sample}" > TMP/jeter.sample


regenie \\
  --step 1 \\
  --bgen ${bgen} \\
  --sample TMP/jeter.sample \\
  --phenoFile ${phenoFile} \\
  --covarFile "${covariates}" \\
  --extract '${keep_rs}' \\
  ${args} \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "${params.prefix}step1"

find ./

"""
}

process ANNOT_VCF {
tag "${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
memory "3G"
input:
	path(vcf)
	path(snpeffDir)
output:
	path("annot.vcf.gz"),emit:output
script:
"""
set -o pipefail

bcftools view "${vcf}" |\\
	snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  eff -dataDir "\${PWD}/${snpeffDir}" \\
                 -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${params.snpeff_db}' |\\
	bcftools view -O z -o annot.vcf.gz

"""
}
