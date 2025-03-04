/** https://github.com/FINNGEN/regenie-pipelines */

include {SNPEFF_DOWNLOAD} from '../../modules/snpeff/download'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../subworkflows/vcf2intervals'
include {BCFTOOLS_CONCAT_CONTIGS} from '../../subworkflows/bcftools/concat.contigs'
include {k1_signature} from '../../modules/utils/k1.nf'
//include {GHOSTSCRIPT_MERGE} from '../..//modules/gs/merge'

workflow {
	if(!file(params.vcf).name.contains(".")) throw new IllegalArgumentException("--vcf missing");
	if(!file(params.samplesheet).name.contains(".")) throw new IllegalArgumentException("--samplesheet missing");
	reference = Channel.of(file(params.fasta), file(params.fai), file(params.dict)).collect()
	user_bed = file(params.bed)

	tobed_ch = JVARKIT_VCF_TO_INTERVALS_01(Channel.fromPath(params.vcf), user_bed)

	snpeff_ch = SNPEFF_DOWNLOAD(params.snpeff_db)

	snsheet_ch = DIGEST_SAMPLESHEET( tobed_ch.bed.
			splitCsv(header:false,sep:'\t').
			map{it[3]}.
			first() ,
		file(params.samplesheet)
		)
	
	wch1_ch = WGSELECT(
			reference,
			snpeff_ch.output,
			snsheet_ch.output,
			tobed_ch.bed.splitCsv(header:false,sep:'\t').
				map{[it[0],((it[1] as int)+1),(it[2] as int),file(it[3]),file(it[4])]}.
				filter{it[0].matches("(chr)?[0-9XY]+")}
			)

	wch1_ch = wch1_ch.output.map{[(it[0].startsWith("chr") ? it[0].substring(3) : it[0] ) , it[1], it[2] ]}


	wch2_ch = BCFTOOLS_CONCAT_CONTIGS(
		wch1_ch ,
		user_bed
		)
	pgen_ch = PLINK2_VCF2PGEN(reference, snsheet_ch.output, wch2_ch.output)
	bgen_ch = PLINK2_MERGE_PGEN(snsheet_ch.output, pgen_ch.output.collect())

	if(!params.covariate.contains(".")) {
		pca_ch = RUN_PCA(snsheet_ch.output, wch2_ch.output)
		make_covar_ch = MAKE_COVARIATES(pca_ch.output)
		covar_ch = make_covar_ch.output
		}
	else
		{
		covar_ch = Channel.of(file(params.covariate))
		}

	/* read header of covariates */
	pca_cols_ch = covar_ch.splitCsv(sep:'\t',header:false).
		take(1).
		flatMap{T->T[2 ..< T.size()]}
	
	/** plot pca for each pair of columns */
	PLOT_PCA(covar_ch, snsheet_ch.output, pca_cols_ch.combine(pca_cols_ch).filter{it[0].compareTo(it[1])<0})

	if(!params.step1_loco.contains(".")) {
		step1 = STEP1( 
			bgen_ch.output /* bgen file */,
			covar_ch /* covariates */, 
			snsheet_ch.output /* samples, sex, phenotype */, 
			pca_ch.output /* contains list of SNP to retain for step 1 */
			)
		loco_ch = step1.output
		}
	else {
		loco_ch = file(params.step1_loco)
		}
	
	
	annot2_ch = MAKE_ANNOT_PER_CTG(wch2_ch.output)

	step2_ch = STEP2(
		bgen_ch.output,
		covar_ch,
		snsheet_ch.output, 
		loco_ch ,
		annot2_ch.output
		)


	ch5 = MERGE_AND_PLOT(
		reference,
		step2_ch.output.collect()
		)
	QQPLOT(ch5.output.flatten().filter{it.name.endsWith(".regenie.gz")})
	MAKE_PDF_ARCHIVE(ch5.output.flatten().filter{it.name.endsWith(".pdf")}.collect())
	}



process DIGEST_SAMPLESHEET {
label "process_single"
tag "${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(vcf)
	path(samplesheet)
output:
	path("pedigree.*"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP
set -x
export LC_ALL=C


# samples in VCF
bcftools query -l "${vcf}" | sort | uniq > TMP/in.vcf.samples.txt
test -s TMP/in.vcf.samples.txt


cat << 'EOF' > TMP/jeter.R

SN= read.table("TMP/in.vcf.samples.txt",header=TRUE,sep="\\t", stringsAsFactors=FALSE,col.names=c("sample"))


T1 <- read.table("${samplesheet}",header=TRUE,sep="\\t", stringsAsFactors=FALSE)

required_columns <- c("sample", "sex", "status")
missing_columns <- setdiff(required_columns, colnames(T1))

if (length(missing_columns) > 0) {
  stop(paste("Erreur : Les colonnes suivantes sont manquantes dans le fichier :", paste(missing_columns, collapse=", ")))
}


if(!"population" %in% colnames(T1)) {
  # Copy the content of 'status' column into the new 'population' column
  T1\$population <- T1\$status
}

T1\$population[T1\$population == ""] <- "NA"
T1\$population[T1\$population == "."] <- "NA"



T1 <- T1[T1\$sample %in% SN\$sample,]

if (any(duplicated(T1\$sample))) {
  stop("DUPLICATE SAMPLES")
}

T2 <- T1[T1\$status=='case',]\$sample
write.table(T2,"TMP/cases.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
T2 <- T1[T1\$status=='control',]\$sample
write.table(T2,"TMP/controls.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)


# Conversion des valeurs de sex pour PLINK (1 = male, 2 = female, 0 = unknown)
T1\$sex <- ifelse(T1\$sex == "male", 1, ifelse(T1\$sex == "female", 2, 0))
T1\$status <- ifelse(T1\$status == "case" , "1", ifelse(T1\$status == "control" , "0" , "NA"))

# header
T1\$FID <- T1\$sample
T1\$IID <- T1\$sample


T2<-T1[,c("sample","population")]
write.table(T2,file="TMP/sample2pop.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

T1 <- T1[, c("FID", "IID", "sex", "${params.status}")]
write.table(T1,file="TMP/plink.ped", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
EOF

R --no-save --vanilla < TMP/jeter.R


sort TMP/cases.txt > pedigree.cases.txt
sort TMP/controls.txt > pedigree.controls.txt
cat TMP/cases.txt TMP/controls.txt > pedigree.all.samples.txt
mv TMP/plink.ped pedigree.plink.ped
mv TMP/sample2pop.txt pedigree.sample2population.tsv
"""
}

process WGSELECT {
label "process_quick"
tag "${contig}:${start1}-${end} ${vcf.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(genome)
	path(snpeffDir)
	path(samplesheet_files)
	tuple val(contig),val(start1),val(end),path(vcf),path(vcf_idx)
output:
	tuple val(contig),path("*.bcf"),path("*.csi"),emit:output
script:
	def interval  = "${contig}:${start1}-${end}"
	def fasta = genome.find{it.name.endsWith("a")}
	def fai = genome.find{it.name.endsWith(".fai")}
	def max_alleles = 3
	def polyx = 10
	def p_assoc = 1E-6
	def f_missing = 0.01
	def gnomad_population = params.gnomad_population?:"AF_nfe"
	def minDP=15
	def maxDP=300
	def all_samples = samplesheet_files.find{it.name.endsWith("all.samples.txt")}
	def controls_file = samplesheet_files.find{it.name.endsWith(".controls.txt")}
	def cases_file = samplesheet_files.find{it.name.endsWith(".cases.txt")}
"""
hostname 1>&2
mkdir -p TMP
set -x
export LC_ALL=C

	# select samples in region
	bcftools view --threads ${task.cpus} \\
		--regions "${interval}" \\
		--samples-file ${all_samples}  \\
		--trim-unseen-allele \\
		--trim-alt-alleles \\
		--no-update \\
		-O u -o TMP/jeter2.bcf "${vcf}"
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	## normalize #################################################################
	bcftools norm --threads ${task.cpus}  -f ${fasta} --multiallelics -any -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	# remove stuff #################################################################
	bcftools annotate -x 'ID,INFO,FILTER,QUAL'  --threads ${task.cpus} -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf 

	# update AC,AN,AF ############################################################
        bcftools +fill-tags --threads ${task.cpus}  -O u -o TMP/jeter2.bcf TMP/jeter1.bcf  -- -t AF,AN,AC
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	## too many alleles
	bcftools view --threads ${task.cpus}  -m 2 -M ${max_alleles} -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	# not in pedigree ###########################################################################
	bcftools view --threads ${task.cpus}  -i 'AC[*]>0' -O u -o TMP/jeter2.bcf TMP/jeter1.bcf 
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	# ignore spanning deletions #################################################################
	bcftools view --threads ${task.cpus} -e 'ALT=="*"' -O u -o TMP/jeter2.bcf TMP/jeter1.bcf 
	mv TMP/jeter2.bcf TMP/jeter1.bcf


	## F_MISSING ##################################################################################
	bcftools view --exclude-uncalled    ${contig.matches("(chr)?Y")?"":"-i \"F_MISSING < ${f_missing}\""}  -O u  -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	## high and low depth ########################################################################
	bcftools view TMP/jeter1.bcf |\\
		jvarkit -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  vcffilterjdk --nocode -e 'final double dp= variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasDP()).mapToInt(G->G.getDP()).average().orElse(${minDP}); if( dp<${minDP} || dp>${maxDP}) return false; return true;' |\\
		bcftools view -O u -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	## polyx ###################################################################################
	if ${ (polyx as int) > 1 } ; then
	
		bcftools view TMP/jeter1.bcf |\\
		jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  vcfpolyx -R "${fasta}" --tag POLYX -n "${polyx}" |\\
		bcftools view -e 'FILTER~"POLYX_ge_${polyx}"' -O u -o TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter1.bcf
	fi

	# bcftools contrast ########################################################################
	bcftools +contrast -0 ${controls_file} -1 ${cases_file} -O u TMP/jeter1.bcf |\\
 		bcftools view -e 'INFO/PASSOC < ${p_assoc}' -O b -o TMP/jeter2.bcf
	mv -v TMP/jeter2.bcf TMP/jeter1.bcf

	# SNPEFF ###################################################################################
        bcftools view TMP/jeter1.bcf |\\
	        snpEff -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  eff -dataDir "\${PWD}/${snpeffDir}" \\
                 -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${params.snpeff_db}' |\\
        	bcftools view -O u -o TMP/jeter2.bcf
        mv TMP/jeter2.bcf TMP/jeter1.bcf

	# CADD #####################################################################################
        if ${params.cadd!=null && !params.cadd.contains(".")} && test -f '${params.cadd}'
        then
                bcftools view TMP/jeter1.bcf |\\
                        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfcadd --tabix '${params.cadd}' | \\
			bcftools view -O u -o TMP/jeter2.bcf
                mv -v TMP/jeter2.bcf TMP/jeter1.bcf
        fi

	# GNOMAD ####################################################################################
        if ${params.gnomad!=null && !params.gnomad.contains(".")} && test -f "${params.gnomad}"
        then
                bcftools view TMP/jeter1.bcf |\\
                        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfgnomad \\
                        --bufferSize 10000 \\
                        --gnomad "${params.gnomad}" \\
                        --fields "${gnomad_population}"  \\
                        --max-af '1.0' |\\
                bcftools view --apply-filters '.,PASS' -O b -o TMP/jeter2.bcf
                mv TMP/jeter2.bcf TMP/jeter1.bcf
        fi


	# a last norm to remove any dup ? ########################################################################
	bcftools norm --threads ${task.cpus}  -f ${fasta} --rm-dup all -O b -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	# rename chromosomes  ####################################################################################
	awk -F '\t' '{printf("%s\t%s\\n",\$1,\$1);}' '${fai}' | sed 's/\tchr/\t/' > TMP/rename.tsv
	bcftools annotate --rename-chrs TMP/rename.tsv --threads ${task.cpus} -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf

	# update variant ID #######################################################################################
	bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' --threads ${task.cpus} -O u -o TMP/jeter2.bcf TMP/jeter1.bcf
	mv TMP/jeter2.bcf TMP/jeter1.bcf


bcftools view --threads ${task.cpus} -O b -o "wgselect.${contig}_${start1}_${end}.bcf" TMP/jeter1.bcf
bcftools index --threads ${task.cpus} "wgselect.${contig}_${start1}_${end}.bcf"
"""
}



process PLINK2_VCF2PGEN {
label "process_short"
afterScript "rm -rf TMP"
tag "chr${contig}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(genome) //for PAR region
	path(pedigree_files)
	tuple val(contig),path(vcf_files)
output:
	path("${contig}.*"),emit:output
script:
	def k1 = k1_signature()
	def fai = genome.find{it.name.endsWith(".fai")}
	def vcf = vcf_files.find{it.name.endsWith(".vcf.gz") || it.name.endsWith(".bcf")}
	def ped = pedigree_files.find{it.name.endsWith(".plink.ped")}
	def args = contig.matches("(chr)?[XY]")?" --update-sex ${ped} --split-par \${HG}":""
"""
set -o pipefail
mkdir -p TMP

cat << EOF | sort -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thg19
1:${k1.hg38}\thg38
chr1:${k1.hg19}\thg19
chr1:${k1.hg38}\thg38
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sort -t '\t' -k1,1 > TMP/jeter.b

HG=`join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b | head -n1`

plink2 ${vcf.name.endsWith(".bcf")?"--bcf":"--vcf"} "${vcf}"  \\
	${args} \\
	--make-pgen erase-phase \\
	--threads ${task.cpus} \\
	--out "${contig}"
"""
}



process PLINK2_MERGE_PGEN {
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(pedigree_files)
        path("INPUT/*")
output:
        path("merged.*"),emit:output
script:
	def ped = pedigree_files.find{it.name.endsWith(".plink.ped")}
"""
set -o pipefail
mkdir -p TMP
set -x
find INPUT -type l -name "*.pgen" | sed 's/\\.pgen\$//' > TMP/jeter.list

if test `wc -l < TMP/jeter.list` -eq 1
then
	ln -s INPUT/*.log merged.log
	ln -s INPUT/*.pvar merged.pvar
	ln -s INPUT/*.psam merged.psam
	ln -s INPUT/*.pgen merged.pgen
else

plink2 --pmerge-list TMP/jeter.list \\
        --make-pgen \\
	--threads ${task.cpus} \\
        --out "merged"
fi


#
# UPDATE file psam with the information of the pedigree
#

# sort on name, keep original order in 1st column
awk -F '\t' '(NR>1) {printf("%d\t%s\\n",NR,\$1);}' merged.psam |\\
	sort -T TMP -t '\t' -k2,2 > TMP/jeter.a

head TMP/jeter.a 1>&2

cut -f2,3,4 "${ped}" |\\
	tail -n +2 |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

head TMP/jeter.b 1>&2


join -t '\t' -1 2 -2 1 -e 'NA' -a 1 -o '1.1,1.2,2.2,2.3' TMP/jeter.a TMP/jeter.b |\\
	sort -t '\t' -k1,1n |\\
	cut -f 2- |\\
	awk -F '\t' 'BEGIN{printf("#FID\tIID\tSEX\t${params.status}\\n");} {printf("%s\t%s\t%s\t%s\\n",\$1,\$1,\$2,\$3);}' > TMP/merged.new.psam

mv  merged.psam original.merged.psam

mv  TMP/merged.new.psam merged.psam
"""
}

process MAKE_COVARIATES {
executor "local"
input:
	path(plink_files)
output:
	path("covariates.tsv"),emit:output
script:
	mds = plink_files.find{it.name.endsWith(".mds")}
"""
awk 'BEGIN{printf("FID\tIID\tY1\tY2\tY3\\n");} (NR>1) {printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$4,\$5,\$6);}' "${mds}" > covariates.tsv
"""
}


workflow TODO {
	bgen = Channel.of(file(params.bgen),file(params.bgen.substring(0,params.bgen.lastIndexOf("."))+".sample")).collect()
	prune_ch = LD_PRUNING(bgen)
	vcf_ch = MAKE_VCF(bgen)

	annot_ch = ANNOT_VCF(vcf_ch.output, snpeff_db.output)

	}

/** STEP 1 of regenie cannot use more than 'X' alleles,  LD PRUNING is required */

process LD_PRUNING {
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	//path(genome)
	path(bgen_files)
output:
	//path("keep.txt"),emit:output
script:
	//def fai = genome.find{it.name.endsWith(".fai")}
	def bgen = bgen_files.find{it.name.endsWith(".bgen")}
	def sample = bgen_files.find{it.name.endsWith(".sample")}

"""
mkdir -p TMP


#### --chr-set `awk -F '\t' '(\$1 ~ /^(chr)?[0-9]+\$/ )' '${fai}' | wc -l `  no-xy no-mt

set -x

plink2 --bgen ${bgen} ref-first \\
	--sample ${sample} \\
	--allow-no-sex \\
	--rm-dup force-first \\
        --indep-pairwise 50 1 0.2 \\
        --out TMP/plink


find TMP 1>&2

plink2 --bgen ${bgen} ref-first \\
	--sample TMP/jeter.sample \\
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
	bcftools annotate -x 'INFO,QUAL,FILTER' -O z -o TMP/sites.vcf.gz

mv TMP/sites.vcf.gz ./

"""
}

process STEP1 {
label "process_quick_high"
conda "${moduleDir}/../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
	path(bgen_files)
        path(covariates)
	path(pheno_files)
	path(plink_files)
output:
	path("${params.prefix}step1_pred.list"),emit:output
	path("${params.prefix}step1.log"),emit:log
script:
	def pgen = bgen_files.find{it.name.endsWith(".pgen")}
	def keep_rs = plink_files.find{it.name.endsWith("keep.id.txt")}
	def args = "--bsize 1000 --bt --phenoCol Y1 "
	def ped = pheno_files.find{it.name.endsWith(".plink.ped")}
"""

mkdir -p TMP/OUT
set -x


regenie \\
  --step 1 \\
  --pgen \$(basename ${pgen} .pgen) \\
  --phenoFile ${ped} \\
  --phenoColList `head -n 1 ${ped} | cut -f4- |tr "\t" ","` \\
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

process MAKE_ANNOT_PER_CTG {
tag "chr${contig}"
label "process_quick_high"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(contig),path(vcf_files)
output:
	tuple val(contig),path("*.annot.txt"),path("*.mask.txt"),path("*.setfile.txt"),emit:output
script:
	def vcf = vcf_files.find{it.name.endsWith(".bcf") || it.name.endsWith(".vcf.gz")}
"""
set -o pipefail
mkdir -p TMP

JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" | head -n1`

cp -v "${moduleDir}/Minikit.java" TMP/
javac -sourcepath TMP -cp "\${JVARKIT_JAR}" -d TMP TMP/Minikit.java

bcftools view -G -O v '${vcf}' |\\
	java -Xmx${task.memory.giga}g \
		-Djava.io.tmpdir=TMP  \\
		-cp TMP:\${JVARKIT_JAR} Minikit \\
			-a TMP/jeter.annot.txt \\
			-s TMP/jeter.setfile.txt \\
			-m TMP/jeter.mask.txt \\
			-w -1


mv -v TMP/jeter.annot.txt "regenie.${contig}.annot.txt"
mv -v TMP/jeter.mask.txt "regenie.${contig}.mask.txt"
mv -v TMP/jeter.setfile.txt "regenie.${contig}.setfile.txt"
"""
}


process STEP2 {
label "process_quick_high"
tag "chr${contig}"
conda "${moduleDir}/../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
        path(bgen_files)
        path(covariates)
        path(pheno_files)
	path(pred_list)
	tuple val(contig),path(annot),path(mask),path(setfile)
output:
        path("*.regenie.gz"),emit:output
script:
	def pgen = bgen_files.find{it.name.endsWith(".pgen")}
	def ped = pheno_files.find{it.name.endsWith(".plink.ped")}
"""

mkdir -p TMP/OUT
set -x

regenie \\
  --step 2 \\
  --pgen \$(basename ${pgen} .pgen) \\
  --phenoFile ${ped} \\
  --covarFile "${covariates}" \\
  --pred ${pred_list} \\
  --mask-def ${mask} \\
  --phenoCol ${params.status} \\
  --bt \\
  --bsize 1000 \\
  --set-list ${setfile} \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "step2.${contig}" \\
  --anno-file ${annot} \\
  --aaf-bins ${params.freq} \\
  --vc-maxAAF ${params.vc_maxAAF} \\
  --bsize 200 \\
  --vc-tests "${params.vc_tests}" \\
  --check-burden-files \\
  --weights-col ${params.weight_column?:4}

gzip --best *.regenie
"""
}


process MERGE_AND_PLOT {
conda "${moduleDir}/../../conda/bioinfo.01.yml"
//afterScript "rm -rf TMP"
label "process_quick"
input:
	path(genome)
	path("INPUT/*")
output:
	path("*.regenie.*"),emit:output
script:

	def dict = genome.find{it.name.endsWith(".dict")}
	def fai = genome.find{it.name.endsWith(".fai")}
	def fasta = genome.find{it.name.endsWith("a")}

"""
mkdir -p TMP


JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" | head -n1`

cp -v "${moduleDir}/Minikit2.java" TMP/Minikit.java
javac -sourcepath TMP -cp "\${JVARKIT_JAR}" -d TMP TMP/Minikit.java

cut -f1,2 '${fai}' | awk -F '\t' '( \$1 ~ /^(chr)?[0-9XY]+\$/ )' > TMP/jeter.fai

# find one header
cat INPUT/*.regenie.gz |\\
	gunzip -c |\\
	grep  '^CHROM' -m1 > TMP/jeter.txt

test -s TMP/jeter.txt

set -o pipefail


cat INPUT/*.regenie.gz |\\
	gunzip -c |\\
	grep -v "^#" |\\
	grep -v '^CHROM' |\\
	LC_ALL=C sort -T TMP -t ' ' -k1,1V -k2,2n |\\
	uniq >> TMP/jeter.txt

java -cp "\${JVARKIT_JAR}:TMP" Minikit -R "${fasta}" -o TMP < TMP/jeter.txt
rm TMP/jeter.txt

cat << '__EOF__' > TMP/jeter.R
mft <- read.table("TMP/manifest.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE, colClasses=c("character","character","character"),)
head(mft)

fai <- read.table("TMP/jeter.fai", header=FALSE, sep="\t", stringsAsFactors=FALSE, colClasses=c("character","double"))
colnames(fai) <- c("contig", "size")

pos2index <- function(CHROM, POS) {  
  total_size <- 0.0
  for (i in 1:nrow(fai)) {
    if (fai\$contig[i] == CHROM) {
      return(total_size + POS)
    }
    total_size <- total_size + fai\$size[i]
  }
    stop("Chromosome non trouve")
}



chrom_positions <- c(0, cumsum(fai\$size))
chrom_colors <- rep(c("azure", "azure2"), length.out=nrow(fai))


last_chrom <- fai\$contig[nrow(fai)]
last_size <- fai\$size[nrow(fai)]
max_genome_size <- pos2index(last_chrom, last_size)


for(i in seq_len(nrow(mft))) {

data <- read.table(mft[i,]\$filename, header=TRUE, sep=" ", stringsAsFactors=FALSE)
data\$x <- mapply(pos2index, data\$CHROM, data\$GENPOS)

fileout <- paste("TMP/freq_",mft[i,]\$freq,"_",mft[i,]\$test_name,"_",mft[i,]\$annot,".regenie.pdf",sep="")

pdf(fileout, width=20, height=6)

x_lim <-  c(0,max_genome_size)
y_lim <- c(0,1+max(data\$LOG10P))

plot(	NULL,
	main=paste("Test:",mft[i,]\$test_name," AF:",mft[i,]\$freq," Annot:",mft[i,]\$annot,sep=""),
	sub= "${params.prefix}regenie",
	xlab="${fasta}",
	ylab="-log10(PVALUE)",
	xlim = x_lim,
	ylim = y_lim
	)

for (i in seq_along(fai\$contig)) {
  rect(chrom_positions[i], par("usr")[3], chrom_positions[i+1], par("usr")[4], col=chrom_colors[i], border=NA)
}

# real plot
points(data\$x,
	data\$LOG10P,
	pch=20,
	col="black",
	xlim = x_lim,
	ylim = y_lim
	)

# vertical bars
abline(v=chrom_positions, col="black", lty=2)
abline(h=6,col="red",lty=2)
abline(h=5,col="green",lty=2)

dev.off()
}
__EOF__

R --no-save < TMP/jeter.R
find TMP -type f

find TMP -type f -name "*.regenie" | while read F
do
	tr " " "\t" <  "\${F}" | bgzip > "\${F}.gz"
	tabix -S 1 -s 1 -b 2 -e 2 "\${F}.gz"
done

mv TMP/*.pdf ./
mv TMP/*.gz ./
mv TMP/*.tbi ./
"""
}

/*
#bgzip TMP/jeter.tsv
#tabix -S 1 -s 1 -b 1 -e 1 TMP/jeter.tsv.gz
#mv TMP/jeter.tsv.gz     regenie.${test_name}.${mask_name}.${status_name}.${freq}.tsv.gz
#mv TMP/jeter.pdf     regenie.${test_name}.${mask_name}.${status_name}.${freq}.pdf
#mv TMP/jeter.tsv.gz.tbi regenie.${test_name}.${mask_name}.${status_name}.${freq}.tsv.gz.tbi
*/


process QQPLOT {
label "process_quick"
tag "${regenie.name}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(regenie)
output:
	path("*.png"),emit:output
script:
"""
mkdir -p TMP
cat << '__EOF__' > TMP/jeter.R
T1 <- read.table(file="${regenie}",sep="\t",header=TRUE, stringsAsFactors=FALSE)

png("${regenie.name}.png")
pvector <- T1\$LOG10P
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector) ))
plot(e, o, pch=20, 
	main="${regenie.name}",
	xlab=expression(Expected~~-log[10](italic(p))), 
	ylab=expression(Observed~~-log[10](italic(p))),
	)
 abline(0,1,col="red")
dev.off()
__EOF__

R --no-save < TMP/jeter.R
"""
}


workflow RUN_PCA {
	take:
		pedigree_files
		rows //[contig,(vcf+vcf_idx)]
	main:
		ch1 = PCA_PER_CONTIG(rows.filter{it[0].matches("(chr)?[0-9]+")})
		all_plink_ch = ch1.output.collect().flatten().collect()
		ch2 = MERGE_PIHAT(all_plink_ch)
		PLOT_AVERAGE_PIHAT(pedigree_files, ch2.output)
	emit:
		output = ch2.output
}


process PCA_PER_CONTIG {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
tag "${contig}"
afterScript "rm -rf TMP"
input:
        tuple val(contig),path(vcf_files)
output:
	path("plink.${contig}.*"),emit:output
script:
	def vcf = vcf_files.find{it.name.endsWith(".bcf") || it.name.endsWith(".vcf.gz")}
	def pihatmaf = 0.1
	def pihatMinGQ = 20
	def f_missing= 0.05
	def minDP= 10
	def maxDP= 300
"""
hostname 1>&2
set -o pipefail
set -x

mkdir -p TMP

echo '${contig}' | sed 's/^chr//' > TMP/chroms.txt

bcftools view -m2 -M2  --apply-filters '.,PASS' --regions `cat TMP/chroms.txt` --types snps -O u "${vcf}" |\
	bcftools view --exclude-uncalled  --min-af "${pihatmaf}" --max-af "${1.0 - (pihatmaf as Double)}"  -i 'AC>0 ${contig.matches("(chr)?Y")?"":"&& F_MISSING < ${f_missing}"}'  -O v |\\
	jvarkit -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  vcffilterjdk --nocode -e 'final double dp= variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasDP()).mapToInt(G->G.getDP()).average().orElse(${minDP}); if( dp<${minDP} || dp>${maxDP}) return false; if (variant.getGenotypes().stream().filter(G->G.isCalled() && G.hasGQ()).anyMatch(G->G.getGQ()< ${pihatMinGQ} )) return false; return true;' |\\
	bcftools annotate  -x '^INFO/AC,INFO/AF,INFO/AN,QUAL,^FORMAT/GT' -O z -o TMP/jeter2.vcf.gz

mv -v TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz
bcftools query -f "\\n" TMP/jeter1.vcf.gz | wc -l 1>&2

# convert VCF to plink (BCF marche pas ?)
plink \\
	--vcf TMP/jeter1.vcf.gz \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--threads ${task.cpus} \\
	--make-bed \\
	--out TMP/jeter1 1>&2

find TMP -type f 1>&2

#
# Dans le fichier hardyweinberg.hwe enlever les variants dont la colonne P < 0.00001
#
##	--allow-no-sex \\

plink \\
	--bfile TMP/jeter1 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--hardy gz \\
	--threads ${task.cpus} \\
	--out "TMP/hardyweinberg.txt" 1>&2

find TMP -type f 1>&2
file TMP/hardyweinberg.txt.hwe.gz 1>&2

gunzip -c TMP/hardyweinberg.txt.hwe.gz |\\
	awk '(\$9 <  0.00001) {print \$2}'  > TMP/xclude_ids.txt

wc -l TMP/xclude_ids.txt 1>&2


# remove variants with those IDS
plink \\
	--bfile TMP/jeter1 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--threads ${task.cpus} \\
	--make-bed \\
	--exclude TMP/xclude_ids.txt \\
	--out TMP/jeter2 1>&2

# rm -v TMP/jeter1*
find TMP -type f 1>&2

#
# selection de SNP independants 
#
plink \\
	--bfile  TMP/jeter2 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--indep-pairwise 50 10 0.2 \\
	--out TMP/jeter3 1>&2

find TMP -type f 1>&2
wc TMP/jeter3.prune.in 1>&2
head TMP/jeter3.prune.in 1>&2


plink \\
	--bfile  TMP/jeter2 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--extract TMP/jeter3.prune.in \\
	--make-bed \\
	--out TMP/jeter4 1>&2

# rm -v TMP/jeter2*
find TMP -type f 1>&2

plink \\
	--bfile  TMP/jeter4 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--r2 --ld-window 50 --ld-window-kb 5000 --ld-window-r2 0.2 \\
	--out TMP/jeter5 1>&2

test -f TMP/jeter5.ld

find TMP -type f 1>&2
awk '{print \$3;}'  TMP/jeter5.ld  | uniq  | sort -T TMP | uniq > TMP/snpInLD.txt

plink \\
	--bfile TMP/jeter4 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--exclude TMP/snpInLD.txt \\
	--make-bed \\
	--out TMP/plink.${contig}.indepSNP_data

#
# it is not recommened to use more than 1000000 variants in step 1 of regenie
#
# find TMP -type f 1>&2
# awk '{print \$2}' TMP/plink.${contig}.indepSNP_data.bim | sort -T TMP | uniq > TMP/${contig}.keep.id.txt

mv -v TMP/plink.${contig}.* ./
"""
}



process MERGE_PIHAT {
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path("PLINK/*")
output:
	path("plink*"),emit:output
script:
	def num_components=3
"""
hostname 2>&1
set -x
mkdir -p TMP
find PLINK/ 1>&2

find PLINK/  -name "*.bim" | sed 's/\\.bim\$//' | sort > TMP/jeter.list
test -s TMP/jeter.list

plink \\
	--merge-list TMP/jeter.list \\
        --allow-extra-chr \\
        --allow-no-sex \\
        --make-bed \\
        --out TMP/jeter1 1>&2

# it is not recommened to use more than 1000000 variants in step 1 of regenie
#
awk 'BEGIN {srand(0);} {printf("%f\t%s\\n", rand(), \$2);}' TMP/jeter1.bim |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1g |\\
	cut -f 2 |\\
	head -n 1000000 > TMP/plink.keep.id.txt

test -s TMP/plink.keep.id.txt

plink \\
	--bfile TMP/jeter1 \\
	--genome \\
	--out TMP/plink 1>&2

find TMP -type f 1>&2

plink \\
	--bfile TMP/jeter1 \\
	--read-genome TMP/plink.genome \\
	--mds-plot ${num_components} \\
	--cluster \\
	-out TMP/plink 1>&2

find TMP -type f 1>&2




mv -v TMP/plink* ./
"""
}

process PLOT_PCA {
tag "${col1} vs ${col2}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(covar)
	path(plink_files)
	tuple val(col1),val(col2)
output:
	path("*.png"),emit:output
script:
	def sample2pop = plink_files.find{it.name.endsWith(".sample2population.tsv")}
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.R
T1 <- read.table(file="${covar}",sep="\\t",header=TRUE, stringsAsFactors=FALSE)

T2 <- read.table(file="${sample2pop}",sep="\\t",header=FALSE, stringsAsFactors=FALSE)
colnames(T2) <- c("IID", "collection")
T1 <- merge(T1, T2, by="IID", all.x=TRUE)
collection_colors <- rainbow(length(unique(T1\$collection)))
names(collection_colors) <- unique(T1\$collection)
T1\$color <- collection_colors[T1\$collection]

png("pihat.${col1}_${col2}.png")
plot(
	x=T1\$${col1},
	y=T1\$${col2},
	col=T1\$color,
	xlab="${col1}",
	ylab="${col2}",
	pch = 19,
	main="PCA: ${col1} x ${col2}",
	sub = "${covar.name}"
	)
legend("topright", legend=names(collection_colors), col=collection_colors, pch=19, title="Collections")
dev.off()
__EOF__

R --no-save < TMP/jeter.R


"""
}


process PLOT_AVERAGE_PIHAT {
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(pedigree_files)
	path(plink_files)
output:
	path("*.png"),emit:output
script:
	def plink_genome = plink_files.find{it.name.endsWith(".genome")}
	def maxPiHat = 0.1
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.R

data <- read.table("${plink_genome}", header = TRUE, stringsAsFactors=FALSE)

# Calculate the average PI_HAT per IID1
unique_IID1 <- unique(data\$IID1)
average_PI_HAT <- sapply(unique_IID1, function(iid1) {
  mean(data\$PI_HAT[data\$IID1 == iid1], na.rm = TRUE)
})

# Create a data frame for the results
T1 <- data.frame(S = unique_IID1, X = average_PI_HAT)

png("sample2avg.pihat.png")
boxplot(
	T1\$X ,
	ylim=c(0,max(T1\$X)),
	main="AVG(PIHAT)/SAMPLE",
	xlab="Sample",
	ylab="pihat"
	)
abline(h=${maxPiHat},col="blue");
dev.off()
__EOF__
R --no-save < TMP/jeter.R
"""
}

process PIHAT_PER_SAMPLE {
tag "${plink_genome}"
label "process_quick"
afterScript "rm -rf TMP"
input:
	path(pedigree_files)
        path(plink_genome)
output:
        path("*.png"),emit:output
script:
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.R
png("${prefix}plot.pihat.png")
genome <- read.table(file="TMP/jeter.tsv",sep="\\t",header=FALSE)
plot(genome\$V2,ylim=c(0,1.0),xlab="Individuals Pair", ylab="PI-HAT", main="PI-HAT")
abline(h=${maxPiHat},col="blue");
dev.off()
__EOF__

R --no-save < TMP/jeter.R
"""
}

process ZOB_TODO {
script:
"""
##
## create MULTIQC CONFIG
##
cat << EOF > "${prefix}multiqc.config.yaml"
custom_data:
  pihat_plot1_${params.step_id}:
    parent_id: pihat_section_${params.step_id}
    parent_name: "PIHAT ${params.step_name}"
    parent_description: "${whatispihat}"
    section_name: "Pihat  ${params.step_name}"
    description: "plot of pihat"
  pihat_plot2_${params.step_id}:
    parent_id: pihat_section_${params.step_id}
    parent_name: "PIHAT  ${params.step_name}"
    parent_description: "${whatispihat}"
    section_name: "Sample to Average pihat"
    description: "Sample to Average pihat"
sp:
  pihat_plot1_${params.step_id}:
    fn: "${prefix}plot.pihat.png"
  pihat_plot2_${params.step_id}:
    fn: "${prefix}plot.sample2avg.pihat.png"
ignore_images: false
EOF



cat << EOF > TMP/${prefix}multiqc.genome_mqc.html
<!--
id: 'highpihat_${params.step_id}'
parent_id: pihat_section_${params.step_id}
parent_name: "PIHAT ${params.step_name}"
parent_description: "${whatispihat}"
section_name: 'High Pihat ${params.step_name}'
description: ' 10 higher pihat pairs.'
-->
EOF

awk '(NR>1) {printf("%s\t%s\t%s\\n",\$1,\$3,\$10);}' TMP/plink.genome |\
	LC_ALL=C sort -t '\t' -T TMP -k3,3gr |\
	head -n 10 |\
	awk -F '\t' 'function sn(SN) {nuscore=split(SN,a,/_/); nuscore = nuscore/4; sample="";for(i=1;i<=nuscore;i++) sample=sprintf("%s%s%s",sample,(i==1?"":"_"),a[i]); return sample;} BEGIN {printf("<table class=\\"table\\"><tr><th>SN1</th><th>SN2</th><th>PIHAT</th></tr>\\n");} {printf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>\\n",sn(\$1),sn(\$2),\$3);} END {printf("</table>\\n");}' >> TMP/${prefix}multiqc.genome_mqc.html


cat << EOF > TMP/${prefix}multiqc.removed_samples_mqc.html
<!--
id: 'rmsamples_${params.step_id}'
parent_id: pihat_section_${params.step_id}
parent_name: "PIHAT ${params.step_name}"
parent_description: "${whatispihat}"
section_name: 'Removed Samples ${params.step_name}'
description: 'Samples that would be removed with max.pihat=${params.pihat.pihat_max}.'
-->
EOF

awk -F '\t' 'BEGIN {printf("<table class=\\"table\\"><tr><th>Sample</th><th>Pihat</th></tr>\\n");} (\$3=="DISCARD") {printf("<tr><td>%s</td><td>%s</td></tr>\\n",\$1,\$2);} END {printf("</table>\\n");}' TMP/samples.keep.status  >> TMP/${prefix}multiqc.removed_samples_mqc.html


gzip --best TMP/plink.genome
mv TMP/plink.genome.gz "${prefix}plink.genome.txt.gz"

mv TMP/${prefix}multiqc.genome_mqc.html ./
mv TMP/${prefix}multiqc.removed_samples_mqc.html ./

"""
}

process MAKE_PDF_ARCHIVE {
label "process_quick"
input:
	path("PDFS/*")
output:
	path("${params.prefix?:""}archive.zip"),emit:output
script:
	def dir="${params.prefix?:""}archive"
"""
hostname 1>&2
mkdir -p "${dir}"

cat << __EOF__ >  "${dir}/index.html" 
<html>
<head>
<meta charset="UTF-8"/>
<title>${params.prefix?:""}Regenie</title>
<script>
var index=0;
var pdfs=[
__EOF__

find PDFS -name "*.pdf" -printf "\\"%f\\"\\n" | LC_ALL=C sort -T . -V | paste -sd ',' >> "${dir}/index.html"
cp -v PDFS/*.pdf ${dir}/

cat << __EOF__ >> "${dir}/index.html"
];

function change(dx) {
	index+=dx;
	index = index%pdfs.length;
	var E = document.getElementById("theimg");
	E.setAttribute("src",pdfs[index]);
	document.getElementById("x").textContent = ""+(index+1)+"/"+pdfs.length;
	}
</script>
</head>
<body>
<button onclick="change(-1)">Prev</button>
<button onclick="change( 1)">Next</button>
<span id="x">1/x</span>
<br/>
<embed id="theimg" src="" width="1000" height="700" />
</body>
</html>
__EOF__

zip -r0  "${params.prefix?:""}archive.zip" "${dir}"
rm -rf '${dir}'
"""
}
