/** https://github.com/FINNGEN/regenie-pipelines */

include {SNPEFF_DOWNLOAD} from '../../modules/snpeff/download'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../subworkflows/vcf2intervals'
include {BCFTOOLS_CONCAT_CONTIGS} from '../../subworkflows/bcftools/concat.contigs'
include {k1_signature} from '../../modules/utils/k1.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'

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
	
	scores_ch = FUNCTIONAL_ANNOTATION_SCORES()
	the_annot_ch = MAKE_FUNCTIONAL_ANNOT_PER_CTG(scores_ch.output, wch2_ch.output).output
	

	if(!params.sliding_windows.isEmpty()) {
		windows_ch = Channel.of("${params.sliding_windows}").
			flatMap{T->{
			def L=[];
			def a = T.split("[,]");
			if(a.size()%2!=0) throw new IllegalArgumentException("not a pair ${params.sliding}");
			for(int i=0;i<a.size();i+=2) {
				L.add([ (a[i+0] as int) , (a[i+1] as int) ]);
				}
			return L;
			}}

		slide_ch = MAKE_SLIDING( wch2_ch.output.combine(windows_ch) )
		the_annot_ch = the_annot_ch.mix(slide_ch.output)
		}


	if(file(params.select_bed).name.contains(".")) {
		make_bed_ch = MAKE_BED(file(params.select_bed),wch2_ch.output)
		the_annot_ch = the_annot_ch.mix(make_bed_ch.output)
		}


	for_step2_ch = the_annot_ch.splitCsv(header:false,sep:'\t',elem:2).
			map{[it[0],it[1],it[2][0],it[2][1],it[2][2],it[2][3]]}



	step2_ch = STEP2(
		bgen_ch.output,
		covar_ch,
		snsheet_ch.output, 
		loco_ch ,
		for_step2_ch
		)




	ch5 = MERGE_AND_PLOT(
		reference,
		step2_ch.output.groupTuple()
		)
	
	ANNOT_HITS(reference,ch5.output.collect())
	MAKE_PNG_ARCHIVE(ch5.images.flatten().collect())
	}

runOnComplete(workflow)

process DIGEST_SAMPLESHEET {
label "queue_quick"
memory "1G"
time "3h"
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
array 100
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
		--apply-filters '.,PASS' \\
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
awk 'BEGIN{printf("FID\tIID\tY1\tY2\tY3\\n");} (NR>1) {printf("%s\t%s\t%s\t%s\t%s\\n",\$2,\$2,\$4,\$5,\$6);}' "${mds}" > covariates.tsv
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
	--const-fid 1 \\
	--allow-no-sex \\
	--rm-dup force-first \\
        --indep-pairwise 50 1 0.2 \\
        --out TMP/plink


find TMP 1>&2

plink2 --bgen ${bgen} ref-first \\
	--sample TMP/jeter.sample \\
	--extract TMP/plink.prune.in \\
	--const-fid 1 \\
	--make-bed \\
	--out TMP/pruned_data

plink --bfile TMP/pruned_data \\
	--const-fid 1 \\
	--allow-extra-chr \\
	--r2 --ld-window 50 \\
	--ld-window-kb 5000 \\
	--ld-window-r2 0.2 \\
	--out TMP/ld


awk '{print \$3;}'  TMP/ld.ld  | uniq  | sort -T . | uniq > TMP/snpInLD.txt


plink --bfile TMP/pruned_data \\
	--const-fid 1 \\
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
	--const-fid 1 \\
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


process FUNCTIONAL_ANNOTATION_SCORES {
executor "local"
output:
	path("scores.tsv"),emit:output
script:
"""
cat << EOF | tr -s " " | tr " " "\t" > scores.tsv
3_prime_UTR_variant     0.1     UTR,UTR3
5_prime_UTR_premature_start_codon_gain_variant  0.2     UTR,UTR5
5_prime_UTR_truncation  0.5     UTR,UTR5
3_prime_UTR_truncation  0.5     UTR,UTR3
5_prime_UTR_variant     0.2     UTR,UTR5
bidirectional_gene_fusion       1.0	.
conservative_inframe_deletion   0.1     protein_altering
conservative_inframe_insertion  0.3     protein_altering
disruptive_inframe_deletion     0.2     protein_altering
disruptive_inframe_insertion    0.2     protein_altering
downstream_gene_variant 0.1     downstream,updownstream
exon_loss_variant       1.0     protein_altering
frameshift_variant      0.4	protein_altering
gene_fusion     0.9     protein_altering
intergenic_region       0.001	.
intragenic_variant      0.01	.
initiator_codon_variant	0.5     protein_altering
intron_variant  0.05    intronic
missense_variant        0.9     protein_altering
non_coding_transcript_exon_variant      0.1     non_coding
non_coding_transcript_variant   0.1     non_coding
splice_acceptor_variant 0.5     protein_altering,splice
splice_donor_variant    0.5     protein_altering,splice
splice_region_variant   0.5     protein_altering,splice
start_retained_variant	0.1      synonymous
start_lost      0.6     protein_altering
stop_gained     0.9     protein_altering
stop_lost       0.6     protein_altering
stop_retained_variant   0.2     synonymous
synonymous_variant      0.1     synonymous
upstream_gene_variant   0.1     upstream,updownstream
EOF
"""

}




process MAKE_FUNCTIONAL_ANNOT_PER_CTG {
tag "chr${contig}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"   
afterScript "rm -rf TMP"
input:
	path(annotations)
        tuple val(contig),path(vcf_files)
output:
        tuple val("functional"),val(contig),path("OUT/manifest.tsv"),emit:output
script:
        def vcf = vcf_files.find{it.name.endsWith(".bcf") || it.name.endsWith(".vcf.gz")}

"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT

bcftools view -O v '${vcf}' |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regeniefunctionalannot \\
		--annotations "${annotations}" \\
		-f ${params.freq} |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regeniemakeannot \\
		-m "${annotations}" \\
		--prefix "chr${contig}chunk" \\
		--reserve 20 \\
		-o \${PWD}/OUT \\
		--gzip \\
		-N 5000
	
"""
}


process MAKE_BED {
tag "chr${contig} ${select_bed.name}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(select_bed)
        tuple val(contig),path(vcf_files)
output:
        tuple val("user_bed"),val(contig),path("OUT/manifest.tsv"),emit:output
script:
        def vcf = vcf_files.find{it.name.endsWith(".bcf") || it.name.endsWith(".vcf.gz")}
	def min_length=300
"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT

bcftools view --regions-file '${select_bed}' -O v '${vcf}' |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regeniebedannot \\
		--bed "${select_bed}" \\
		--min-length ${min_length} \\
		-f ${params.freq} |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regeniemakeannot \\
		--prefix "chr${contig}_bed_chunk" \\
		-o \${PWD}/OUT \\
		--reserve 10 \\
		--gzip \\
		-N 5000
"""
}

process MAKE_SLIDING {
tag "chr${contig} ${win_size}/${win_shift}"
label "process_quick"
conda "${moduleDir}/../../conda/bioinfo.01.yml"   
afterScript "rm -rf TMP"
input:
        tuple val(contig),path(vcf_files),val(win_size),val(win_shift)
output:
        tuple val("sliding_${win_size}_${win_shift}"),val(contig),path("OUT/manifest.tsv"),emit:output
script:
        def vcf = vcf_files.find{it.name.endsWith(".bcf") || it.name.endsWith(".vcf.gz")}

"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT

bcftools view -O v '${vcf}' |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regenieslidingannot \\
		--window-size "${win_size}" \\
		--window-shift "${win_shift}" \\
		-f ${params.freq} |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regeniemakeannot \\
		--prefix "chr${contig}_${win_size}_${win_shift}_chunk" \\
		-o \${PWD}/OUT \\
		--reserve 20 \\
		--gzip \\
		-N 5000
	
"""
}



process STEP2 {
label "process_quick_high"
array 100
tag "chr${contig} ${title} ${annot.name}"
conda "${moduleDir}/../../conda/regenie.yml"
afterScript "rm -rf TMP"
input:
        path(bgen_files)
        path(covariates)
        path(pheno_files)
	path(pred_list)
	tuple val(title),val(contig),path(annot),path(setfile),path(mask),path(aff)
output:
        tuple val(title),path("*.regenie.gz"),emit:output
script:
	def pgen = bgen_files.find{it.name.endsWith(".pgen")}
	def ped = pheno_files.find{it.name.endsWith(".plink.ped")}
"""

mkdir -p TMP/OUT
set -x

gunzip -c "${annot}" > TMP/annot.txt
gunzip -c "${setfile}" > TMP/setfile.txt
gunzip -c "${mask}" > TMP/mask.txt
gunzip -c "${aff}" > TMP/aaf.txt

regenie \\
  --step 2 \\
  --pgen \$(basename ${pgen} .pgen) \\
  --phenoFile ${ped} \\
  --covarFile "${covariates}" \\
  --pred ${pred_list} \\
  --mask-def TMP/mask.txt \\
  --set-list TMP/setfile.txt \\
  --anno-file TMP/annot.txt \\
  --aaf-file TMP/aaf.txt \\
  --phenoCol ${params.status} \\
  --bt \\
  --bsize 1000 \\
  --lowmem \\
  --lowmem-prefix TMP/regenie_tmp_preds \\
  --threads ${task.cpus} \\
  --out "step2.${contig}" \\
  --aaf-bins ${params.freq} \\
  --vc-maxAAF ${params.vc_maxAAF} \\
  --bsize 200 \\
  --vc-tests "${params.vc_tests}" \\
  --check-burden-files \\
  --weights-col ${params.weight_column?:4} \\
  --write-mask-snplist \\
  --firth --approx \\
  --pThresh 0.01

gzip --best *.regenie
"""
}


process MERGE_AND_PLOT {
tag "${title}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
label "process_quick"
input:
	path(genome)
	tuple val(title),val(L) //path("INPUT/*")
output:
	path("${title}.results.tsv.gz"),emit:output
	path("*.png"),emit:images
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

# do not use chrY because regenie merge it with chrX (see doc)
cut -f1,2 '${fai}' | awk -F '\t' '( \$1 ~ /^(chr)?[0-9X]+\$/ )' > TMP/jeter.fai


cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

##cat INPUT/*.regenie.gz 

xargs -a TMP/jeter.list -L 1 gunzip -c |\\
	grep  '^CHROM' -m1 > TMP/jeter.txt

test -s TMP/jeter.txt

set -o pipefail

#cat INPUT/*.regenie.gz |\\
#	gunzip -c


xargs -a TMP/jeter.list -L 1 gunzip -c |\\
	grep -v "^#" |\\
	grep -v '^CHROM' |\\
	LC_ALL=C sort -T TMP -t ' ' -k1,1V -k2,2n |\\
	uniq >> TMP/jeter.txt

cat TMP/jeter.txt | tr " " "\t" | gzip --best > TMP/${title}.results.tsv.gz

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
data\$LOG10P <-  as.double(data\$LOG10P)

fileout <- paste("TMP/${title}.freq_",mft[i,]\$freq,"_",mft[i,]\$test_name,"_",mft[i,]\$annot,".regenie.manhattan.png",sep="")

png(fileout, width=1500, height=500,units="px")

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

for (k in seq_along(fai\$contig)) {
  rect(chrom_positions[k], par("usr")[3], chrom_positions[k+1], par("usr")[4], col=chrom_colors[k], border=NA)
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

pvector <- data\$LOG10P
# I don't understand why running pow(10,x) and then log10(y) gives a good chart. I hate R.
pvector <- 10^(-pvector)
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector>0]
if(length(pvector) > 0) {
fileout <- paste("TMP/${title}.freq_",mft[i,]\$freq,"_",mft[i,]\$test_name,"_",mft[i,]\$annot,".regenie.qqplot.png",sep="")
png(fileout, width=500, height=500,units="px")


# Observed and expected
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector) ))

plot(x=e,
	y=o,
	pch=20, 
	xlim=c(0, max(e)),
	ylim=c(0, max(o)), 
        xlab=expression(Expected~~-log[10](italic(p))), 
        ylab=expression(Observed~~-log[10](italic(p))),
	main=paste("Test:",mft[i,]\$test_name," AF:",mft[i,]\$freq," Annot:",mft[i,]\$annot,sep=""),
	sub= "${params.prefix}regenie"
	)

abline(0,1,col="red")
dev.off()
}
}
__EOF__

R --no-save < TMP/jeter.R

#
# NO NEED to saved the regenie.gz, and regenie.gz.tbi as I wrote jvarkit+swing+regenie
#
find TMP -type f
#find TMP -type f -name "*.regenie" | while read F
#do
#	tr " " "\t" <  "\${F}" | bgzip > "TMP/${title}.\${F##*/}.gz"
#	tabix -S 1 -s 1 -b 2 -e 2 "TMP/${title}.\${F##*/}.gz"
#done

mv TMP/*.png ./
mv TMP/*.gz ./
"""
}

/*
#bgzip TMP/jeter.tsv
#tabix -S 1 -s 1 -b 1 -e 1 TMP/jeter.tsv.gz
#mv TMP/jeter.tsv.gz     regenie.${test_name}.${mask_name}.${status_name}.${freq}.tsv.gz
#mv TMP/jeter.pdf     regenie.${test_name}.${mask_name}.${status_name}.${freq}.pdf
#mv TMP/jeter.tsv.gz.tbi regenie.${test_name}.${mask_name}.${status_name}.${freq}.tsv.gz.tbi
*/


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
	--const-fid 1 \\
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
	--const-fid 1 \\
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
	--const-fid 1 \\
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
	--const-fid 1 \\
	--allow-no-sex \\
	--indep-pairwise 50 10 0.2 \\
	--out TMP/jeter3 1>&2

find TMP -type f 1>&2
wc TMP/jeter3.prune.in 1>&2
head TMP/jeter3.prune.in 1>&2


plink \\
	--bfile  TMP/jeter2 \\
	--allow-extra-chr \\
	--const-fid 1 \\
	--allow-no-sex \\
	--extract TMP/jeter3.prune.in \\
	--make-bed \\
	--out TMP/jeter4 1>&2

# rm -v TMP/jeter2*
find TMP -type f 1>&2

plink \\
	--bfile  TMP/jeter4 \\
	--const-fid 1 \\
	--allow-extra-chr \\
	--allow-no-sex \\
	--r2 --ld-window 50 --ld-window-kb 5000 --ld-window-r2 0.2 \\
	--out TMP/jeter5 1>&2

test -f TMP/jeter5.ld

find TMP -type f 1>&2
awk '{print \$3;}'  TMP/jeter5.ld  | uniq  | sort -T TMP | uniq > TMP/snpInLD.txt

plink \\
	--bfile TMP/jeter4 \\
	--const-fid 1 \\
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
	--const-fid 1 \\
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
	path("average_pihat.tsv")
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
# save in a file
write.table(T1,"average_pihat.tsv",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)

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

process MAKE_PNG_ARCHIVE {
label "process_quick"
input:
	path("PNGS/*")
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
var index=-1;
var images=[
__EOF__

find PNGS/ -name "*.manhattan.png" -printf "\\"%f\\"\\n" |\\
	sed 's/.manhattan.png"/"/' |\\
	LC_ALL=C sort -T . -V |\\
	paste -sd ',' >> "${dir}/index.html"
cp -v PNGS/*.png ${dir}/

cat << __EOF__ >> "${dir}/index.html"
];

function goTo(idx) {
	if(index==idx || images.length==0) return;
        if(idx<0) idx=images.length-1;
        if(idx>=images.length) idx=0;
        index=idx;
        var E = document.getElementById("theimg");
        E.setAttribute("src",images[index]+".manhattan.png");
        E = document.getElementById("qq");
        E.setAttribute("src",images[index]+".qqplot.png");
        document.getElementById("x").textContent = ""+(index+1)+"/"+images.length;
	E=document.getElementById("select");
	if(E.selectedIndex!=index) E.selectedIndex=index;
        }

addEventListener("load", (event) => {
        goTo(0);
        var E=  document.getElementById("select");
        E.addEventListener("change",()=>{goTo(parseInt(E.value));});
        for(var i in images) {
                var C = document.createElement("option");
                C.setAttribute("value",i);
                C.appendChild(document.createTextNode(images[i]));
                E.appendChild(C);
                }
        });
</script>
</head>
<body>
<button onclick="goTo(index-1)">Prev</button>
<select id="select"></select>
<button onclick="goTo(index+1)">Next</button>
<span id="x">1/x</span>
<br/>
<img id="theimg" alt="manhattan" src=""/><img id="qq" alt="qqplot" src=""/><br/>
</body>
</html>
__EOF__

zip -r0  "${params.prefix?:""}archive.zip" "${dir}"
rm -rf '${dir}'
"""
}

process ANNOT_HITS {
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
label "process_quick"
input:
	path(genome)
	path("INPUT/*")
output:
	path("digest.tsv"),emit:output
script:
	def fai = genome.find{it.name.endsWith(".fai")}
	def dict = genome.find{it.name.endsWith(".dict")}
	def treshold=(params.digest_treshold as double)
        def k1 = k1_signature()
	def n_try = "--tries 10 --waitretry=120"
"""
export LC_ALL=C
set -o pipefail
set -x
mkdir -p TMP

cat << EOF | sort -t '\t' -k1,1 > TMP/jeter.a
1:${k1.hg19}\thttp://grch37.ensembl.org/biomart/martservice
1:${k1.hg38}\thttp://ensembl.org/biomart/martservice
EOF


awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -t '\t' -k1,1 > TMP/jeter.b

BIOMART=`join -t '\t' -1 1 -2 1 -o 1.2 TMP/jeter.a TMP/jeter.b | head -n1`


gunzip -c INPUT/*.gz |\\
	awk -F '\t' '(\$1=="CHROM" || (\$12!="NA" && \$12>${treshold}))' > TMP/jeter.a

awk  -F '\t' '(\$1=="CHROM")' TMP/jeter.a | head -n 1 > TMP/save.header

awk -F '\t' '(\$1!="CHROM") {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$2)-1,\$2,\$0);}' TMP/jeter.a |\\
	jvarkit bedrenamechr -R "${dict}" --column 1 --convert RETURN_ORIGINAL |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n |\\
		uniq > TMP/results.bed
		

cat << __EOF__ > TMP/biomart.01.xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
		<Attribute name = "chromosome_name" />
		<Attribute name = "start_position" />
		<Attribute name = "end_position" />
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "external_gene_name" />
		<Attribute name = "description" />
		<Attribute name = "definition_1006" />
	</Dataset>
</Query>
__EOF__

wget ${n_try} -O TMP/jeter.biomart.tsv "\${BIOMART}?query=\$(cat TMP/biomart.01.xml)"


awk -F '\t' '{OFS="\t";\$2=int(\$2)-1;print;}' TMP/jeter.biomart.tsv |\\
	sed \\
		-e 's/\\[Source:[A-Za-z_0-9 ;:]*\\]//' \\
		-e 's/A location, relative to cellular compartments and structures, occupied by a macromolecular machine. There are three types of cellular components .*//' \
		-e 's/A biological process is the execution of a genetically-encoded biological .*//' |\\
	uniq |\\
	jvarkit bedrenamechr -R "${dict}" --column 1 --convert SKIP |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n -k3,3n -k4,4 > TMP/jeter.biomart.bed


awk -F '\t' '{printf("%s",\$0);for(i=0;i< 7-NF;i++) printf("\t.");printf("\\n");}' < TMP/jeter.biomart.bed |\\
	datamash  --filler=NA --no-strict  -g 1,2,3,4 -c ','  -t '\t' unique 5 unique 6 unique 7 |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n |\\
	uniq > TMP/annot.bed

head TMP/annot.bed 1>&2

tr "\\n" "\t" < TMP/save.header > TMP/digest.tsv

echo -e "gene_chrom\tgene_start\tgene_end\tgene_id\tgene_name\tgene_desc\tgene_GO\tdistance_to_gene" >> TMP/digest.tsv


bedtools closest \\
	-a TMP/results.bed \\
	-b TMP/annot.bed \\
	-d  |\\
	cut -f 4- |\\
	sort -t '\t' -T TMP -k12,12gr >> TMP/digest.tsv


mv TMP/digest.tsv digest.tsv
"""
}
