
include {runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'


def LOAD_DIFFBIND = """

set +u
source \${HOME}/.bash_profile
micromamba activate DIFFBIND
set -u

"""


workflow {
	DIFFBIND(Channel.fromPath(params.samplesheet))
	}

runOnComplete(workflow);


workflow DIFFBIND {
	take:
		samplesheet
	main:
		rows1_ch = samplesheet.splitCsv(header:true, sep:'\t').
			map{T->(!T.containsKey("sample") && T.containsKey("sampleId") ?T.plus("sample":T.sampleId):T)}.
			map{T->(!T.containsKey("bam") && T.containsKey("bamReads") ?T.plus("bam":T.bamReads):T)}

		rows2_ch = SAMTOOLS_IDXSTATS(rows1_ch)
		rows3_ch = rows2_ch.splitText().
			map{T->T[0].plus("countReads":(T[1] as long))}
		min_count_reads = rows3_ch.map{it.countReads}.min()
		rows4_ch = rows3_ch.combine(min_count_reads).
			map{T->T[0].plus("minCountReads":T[1])}
		
		downsample_bam_ch = DOWNSAMPLE_BAM(rows4_ch)


		dba1 = DBA_LOAD(
			downsample_bam_ch.output.
				map{T->T[0].plus("bam":T[1])}.
				collect()
			)
		
		
		dba_count_ch = DBA_COUNT(dba1.output)


		PLOT_DBA(
			dba1.output.map{[it,"OccupancyAnalysis"]}.
			mix(dba_count_ch.output.map{[it,"ClusteringAffinity"]})
			)


	
		cols_ch = Channel.from('Tissue','Condition','Factor').
			combine(dba1.samplesheet.splitCsv(header:true, sep:'\t')).
			filter{T->T[1].containsKey(T[0])}.
			map{T->[T[0],T[1][T[0]]]}.
			unique().
			groupTuple().
			filter{T->T[1].size()>1}

		PLOT_PCA(dba_count_ch.output, cols_ch)


		norm_ch = DBA_NORM(dba_count_ch.output)
		contrast_ch = DBA_CONTRAST(norm_ch.output, cols_ch)
		analyze_ch = DBA_ANALYZE(contrast_ch.output)
		report_ch = DBA_REPORT(analyze_ch.output)

	}

process SAMTOOLS_IDXSTATS {
tag "${row.sample}"
input:
	val(row)
output:
	tuple val(row),path("count.txt"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

samtools idxstats "${row.bam}" | awk -F '\t' 'BEGIN {N=0;} (\$1!="*") {N+=int(\$3);} END {printf("%d",N);}' > count.txt
"""
}


process DOWNSAMPLE_BAM {
tag "${row.sample}"
cpus 4
input:
	val(row)
output:
	tuple val(row),path("${row.sample}.bam"),path("${row.sample}.bam"),emit:output
script:
	def factor = row.countReads> row.minCountReads ? " -s "+(row.minCountReads/(double)row.countReads) : ""
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

samtools view ${factor} --threads ${task.cpus} --write-index  -O BAM -o "${row.sample}.bam##idx##${row.sample}.bam.bai" "${row.bam}"
"""
}



process DBA_LOAD {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(L)
output:
	path("dba_DBA.RData"),emit:output
	path("samplesheet.tsv"),emit:samplesheet
script:
"""
mkdir -p TMP
set +u
source \${HOME}/.bash_profile
micromamba activate DIFFBIND
set -u


cat << EOF | paste -s -d '\t' > TMP/samplesheet.tsv
sampleId
Peaks
PeakCaller
bamReads
Tissue
Condition
QC
Factor
Replicate
EOF


cat << __EOF__ >> TMP/samplesheet.tsv
${L.collect{[
	it.sample,
	it.Peaks,
	(it.PeakCaller?:"narrow"),
	it.bam,
	(it.Tissue?:"NA"),
	(it.Condition?:"NA"),
	(it.QC?:"NA"),
	(it.Factor?:"NA"),
	(it.Replicate?:"1")
	].join("\t")}.join("\n")}
__EOF__

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)

samplesheet = read.csv("TMP/samplesheet.tsv",h=TRUE ,sep="\t", stringsAsFactor = FALSE)
head(samplesheet)

DBA <- dba(sampleSheet=samplesheet)
dba.save(DBA, file='DBA', dir='.', pre='dba_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R

mv TMP/samplesheet.tsv ./
"""
}

process PLOT_DBA {
tag "${dba.name} ${title}"
afterScript "rm -rf TMP"
input:
	tuple path(dba),val(title)
output:
	path("${title}.png"),emit:output
script:
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
png("${title}.png", width= 840, height = 840)
plot(DBA,
	cexRow = 1.4,
	cexCol = 1.4
	)
dev.off()
__EOF__

R --vanilla < TMP/jeter.R
"""
}

/**

The next step is to calculate a binding matrix with scores based on read counts for every
sample (affinity scores), rather than confidence scores for only those peaks called in a specific
sample (occupancy scores).

*/
process DBA_COUNT {
memory "50G"
input:
	path(dba)
output:
	path("dba_count_DBA.RData"),emit:output
script:
"""
mkdir -p TMP
set +u
source \${HOME}/.bash_profile
micromamba activate DIFFBIND
set -u

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
DBA <- dba.count(DBA , bParallel=FALSE)
dba.save(DBA, file='DBA', dir='.', pre='dba_count_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}



process PLOT_PCA {
tag "${colName} ${L.join("|")}"
afterScript "rm -rf TMP"
input:
	path(dba)
	tuple val(colName),val(L)
output:
	path("PCA_${colName.toUpperCase()}.png"),emit:output
script:
"""
mkdir -p TMP

${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
png("PCA_${colName.toUpperCase()}.png", width= 840, height = 840)
dba.plotPCA(DBA,
	attributes=DBA_${colName.toUpperCase()},
	title="${L.join(" vs ")}"
	)
dev.off()
__EOF__

R --vanilla < TMP/jeter.R
"""
}

/**
The next step is to tell DiffBind how the data are to be normalized.

*/
process DBA_NORM {
tag "${dba.name}"
afterScript "rm -rf TMP"
input:
        path(dba)
output:
        path("dba_norm_DBA.RData"),emit:output
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
DBA <- dba.normalize(DBA)
head(DBA)
dba.save(DBA, file='DBA', dir='.', pre='dba_norm_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}


/**  Sets up contrasts for differential binding affinity analysis  */
process DBA_CONTRAST {
tag "${dba.name} ${colName} ${colVals.join("/")}"
afterScript "rm -rf TMP"
input:
        path(dba)
	tuple val(colName),val(colVals)
output:
        tuple val(colName),val(colVals),path("dba_contrast_${colName}_DBA.RData"),emit:output
	path("show.tsv"),optional:true,emit:show
when:
	colVals.size()>1
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
DBA <- dba.contrast(DBA , categories = DBA_${colName.toUpperCase()})

df <- dba.show(DBA, bContrasts = TRUE)

if(!is.null(df)) {
	head(df)
	write.table(df,file= "show.tsv", quote=F, row.names=F, col.names=T, sep="\t")
	}


dba.save(DBA, file='DBA', dir='.', pre='dba_contrast_${colName}_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}

process DBA_ANALYZE {
tag "${dba.name} ${colName} ${colVals.join("/")}"
afterScript "rm -rf TMP"
input:
	tuple val(colName),val(colVals),path(dba)
output:
        tuple val(colName),val(colVals),path("dba_analyze_${colName}_DBA.RData"),emit:output
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
DBA <- dba.analyze(DBA , method=c(DBA_DESEQ2))

#dba.plotVenn(DBA,contrast=1,method=DBA_DESEQ2)
#dba.plotMA(DBA, method=DBA_DESEQ2)
dba.plotMA(DBA, bXY=TRUE)
dba.plotBox(DBA)



dba.save(DBA, file='DBA', dir='.', pre='dba_analyze_${colName}_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}

process DBA_REPORT {
tag "${dba.name} ${colName} ${colVals.join("/")}"
afterScript "rm -rf TMP"

input:
	tuple val(colName),val(colVals),path(dba)
//output:
//        tuple val(colName),path("dba_analyze_DBA_${colName}_.RData"),emit:output
when:
	colVals.size()>1
script:
	def FDR  = 0.1
	def fold = 1
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
DBA.report <- dba.report(DBA , th = 1, bCalled=TRUE, method=DBA_DESEQ2)
DBA.report <- sort(DBA.report)

head(DBA.report)


  df <- data.frame(CHROM=seqnames(DBA.report),
                   starts=start(DBA.report)-1,
                   ends=end(DBA.report),
                   names=c(paste0("peak",1:length(DBA.report))),
                   scores=c(rep("0", length(DBA.report))),
                   strands=c(rep(".", length(DBA.report))),
                   conc = DBA.report\$Conc,
                   conc_cond1=mcols(DBA.report)[,2],
                   conc_cond2=mcols(DBA.report)[,3],
                   fold = DBA.report\$Fold,
                   pval = DBA.report\$"p-value",
                   FDR = DBA.report\$FDR
  		)

  colnames(df)[8] <- paste0("${colVals[0]}","_conc")
  colnames(df)[9] <- paste0("${colVals[1]}","_conc")
  
  write.table(df,file = "TMP/jeter1.bed", quote=F, col.names=T, row.names = F, sep="\t")
  
  temp = as.character(df\$names)
  temp[which(df\$fold <= 0)] <- paste0("${colVals[0]}","_", temp[which(df\$fold <= 0)])
  temp[which(df\$fold >  0)] <- paste0("${colVals[0]}","_", temp[which(df\$fold > 0)])
  
  df\$names=temp
  ###Selecting based on asb(fold change) > 1 and FDR < 1%
  df_select <- df[which(abs(df\$fold) >= ${fold}),]
  
  df_select <- df_select[which(df_select\$FDR <= ${FDR}),]
  
  write.table(df_select[,c(1,2,3,4)],file= "TMP/jeter2.bed", quote=F, row.names=F, col.names=F, sep="\t")


__EOF__

R --vanilla < TMP/jeter.R

module load htslib
sed 's/^CHROM/#CHROM/' TMP/jeter1.bed |\\
	LC_ALL=C sort -t '\t' -T TMP -k1,1 k2,2n |\\
	bgzip > "TMP/jeter1.bed.gz"
tabix --comment '#' -f -p bed TMP/jeter1.bed.gz

mv TMP/jeter1.bed.gz "${colVals[0]}_${colVals[1]}_full.bed.gz"
mv TMP/jeter1.bed.gz.tbi "${colVals[0]}_${colVals[1]}_full.bed.gz.tbi"



sed 's/^CHROM/#CHROM/' TMP/jeter2.bed |\\
	LC_ALL=C sort -t '\t' -T TMP -k1,1 k2,2n |\\
	bgzip > "TMP/jeter2.bed.gz"
tabix --comment '#' -f -p bed TMP/jeter2.bed.gz

mv TMP/jeter2.bed.gz "${colVals[0]}_${colVals[1]}_fdr-${FDR}_fold-${fold}.bed.gz"
mv TMP/jeter2.bed.gz.tbi "${colVals[0]}_${colVals[1]}_fdr-${FDR}_fold-${fold}.bed.gz.tbi"


"""
}
