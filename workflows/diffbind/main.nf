
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


		if(params.downsample as boolean) {

			rows1_ch.map{T->(!T.containsKey("sample") && T.containsKey("SampleId") ?T.plus("sample":T.SampleId):T)}

			idxstats_ch = SAMTOOLS_IDXSTATS(rows1_ch).output.
				splitText().
                                map{T->T[0].plus("countReads":(T[1] as long))}

			//mergeidx_ch = MERGE_IDXSTATS(idxstats_ch.collect())

			min_count_reads = idxstats_ch.map{it.countReads}.min()

			rows4_ch = idxstats_ch.combine(min_count_reads).
				map{T->T[0].plus("minCountReads":T[1])}
		
			downsample_bam_ch = DOWNSAMPLE_BAM(rows4_ch).output.
					map{T->T[0].plus("bam":T[1])}
			}
		else
			{
			downsample_bam_ch = rows1_ch
			}


		dba1 = DBA_LOAD( downsample_bam_ch.collect() )
		
		
		dba_count_ch = DBA_COUNT(dba1.output)


		plot_dba_ch = PLOT_DBA(
			dba1.output.map{[it,"OccupancyAnalysis"]}.
			mix(dba_count_ch.output.map{[it,"ClusteringAffinity"]})
			)


	
		pca_ch = PLOT_PCA(
			dba_count_ch.output
			)


		norm_ch = DBA_NORM(dba_count_ch.output)
		contrast_ch = DBA_CONTRAST(norm_ch.output)
		analyze_ch = DBA_ANALYZE(contrast_ch.output)
		report_ch = DBA_REPORT(analyze_ch.output)
		

		

		plot_peaks_ch  = PLOT_PEAKS(
			rows1_ch.map{T->[(T.SampleID+" "+(T.ncbiSrs?:"")+" "+(T.ncbiSrx?:"")).trim().replace(' ','_'), T.Condition, T.bigWig].join("\t")}.collect(),
			report_ch.filtered.splitCsv(header:true,sep:'\t')
			)

		occupancy_ch = plot_dba_ch.output.filter{F->F.name.contains("OccupancyAnalysis")}
		affinity_ch = plot_dba_ch.output.filter{F->F.name.contains("ClusteringAffinity")}

		pdf_ch = MERGE_PDFS(plot_peaks_ch.output.collect())

		README(
			dba1.samplesheet,
			occupancy_ch,
			affinity_ch,
			pca_ch.output,
			analyze_ch.ma,
			analyze_ch.boxplot,
			report_ch.output,
			pdf_ch.output
			)

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

samtools idxstats "${row.bam}" | awk -F '\t' 'BEGIN {N=0;} (\$1 ~ /^(chr)?[0-9]+/) {N+=int(\$3);} END {printf("%d",N);}' > count.txt
"""
}

process MERGE_IDXSTATS {
executor "local"
tag "N=${L.size()}"
input:
	val(L)
output:
	path("output.tsv"),emit:output
script:
"""
cat << EOF > output.tsv
${L.collect{T->[T.sample,T.bam,T.countReads].join("\t")}.join("\n")}
EOF
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
	path("${params.prefix}samplesheet.tsv"),emit:samplesheet
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

mv TMP/samplesheet.tsv ./${params.prefix}samplesheet.tsv
"""
}

process PLOT_DBA {
tag "${dba.name} ${title}"
afterScript "rm -rf TMP"
input:
	tuple path(dba),val(title)
output:
	path("${params.prefix}${title}.png"),emit:output
script:
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
png("${params.prefix}${title}.png",width= 840, height = 840)
plot(DBA,
	main="${title}",
	sub="${dba.name}",
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
tag "${dba.name}"
afterScript "rm -rf TMP"
input:
	path(dba)
output:
	path("${params.prefix}PCA_CONDITION.png"),emit:output
script:
"""
mkdir -p TMP

${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
png("${params.prefix}PCA_CONDITION.png", width= 840, height = 840)
dba.plotPCA(DBA,
	attributes=DBA_CONDITION,
	title="Condition"
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
script:
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')
DBA <- dba.normalize(DBA)
dba.save(DBA, file='DBA', dir='.', pre='dba_norm_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}


/**  Sets up contrasts for differential binding affinity analysis  */
process DBA_CONTRAST {
tag "${dba.name}"
afterScript "rm -rf TMP"
input:
        path(dba)
output:
      	path("dba_contrast_DBA.RData"),emit:output
	path("show.tsv"),optional:true,emit:show
script:
	/**  By default, dba.contrast won't create contrasts with less than 3 members in each group.  The easiest solution is to set minMembers=2 when calling dba.contrast */
	def minMembers=2
"""
mkdir -p TMP
${LOAD_DIFFBIND}

cp -v  "${dba}" TMP/dba_DBA.RData

cat << '__EOF__' > TMP/jeter.R
library(DiffBind)
DBA <- dba.load(file='DBA', dir='TMP', pre='dba_', ext='RData')

DBA <- dba.contrast(DBA , categories = DBA_CONDITION, minMembers=${minMembers})

df <- dba.show(DBA, bContrasts = TRUE)

if(!is.null(df)) {
	head(df)
	write.table(df,file= "show.tsv", quote=F, row.names=F, col.names=T, sep="\t")
	}


dba.save(DBA, file='DBA', dir='.', pre='dba_contrast_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}

process DBA_ANALYZE {
tag "${dba.name}"
afterScript "rm -rf TMP"
input:
	path(dba)
output:
        path("dba_analyze_DBA.RData"),emit:output
	path("${params.prefix}ma.png"),emit:ma
	path("${params.prefix}boxplot.png"),emit:boxplot
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

png("${params.prefix}ma.png")
dba.plotMA(DBA, bXY=TRUE)
dev.off()

png("${params.prefix}boxplot.png")
dba.plotBox(DBA)
dev.off()


dba.save(DBA, file='DBA', dir='.', pre='dba_analyze_', ext='RData')
__EOF__

R --vanilla < TMP/jeter.R
"""
}

process DBA_REPORT {
tag "${dba.name}"
afterScript "rm -rf TMP"

input:
	path(dba)
output:
	path("*{bed.gz,bed.gz.tbi}"),emit:output
	path("filtered.tsv"),emit:filtered
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


  df <- data.frame(chrom=seqnames(DBA.report),
                   start=start(DBA.report)-1,
                   end=end(DBA.report),
                   name=c(paste0("peak",1:length(DBA.report))),
                   scores=c(rep("0", length(DBA.report))),
                   strands=c(rep(".", length(DBA.report))),
                   conc = DBA.report\$Conc,
                   conc_condition1=mcols(DBA.report)[,2],
                   conc_condition2=mcols(DBA.report)[,3],
                   fold = DBA.report\$Fold,
                   p_value = DBA.report\$"p-value",
                   FDR = DBA.report\$FDR
  		)

  
  write.table(df,file = "TMP/jeter1.bed", quote=F, col.names=T, row.names = F, sep="\t")
  
  ###Selecting based on asb(fold change) > 1 and FDR < 1%
  df_select <- df[which(abs(df\$fold) >= ${fold}),]
  
  df_select <- df_select[which(df_select\$FDR <= ${FDR}),]
  
  write.table(df_select,file= "TMP/jeter2.bed", quote=F, row.names=F, col.names=T, sep="\t")


__EOF__

R --vanilla < TMP/jeter.R

module load htslib
sed 's/^chrom/#chrom/' TMP/jeter1.bed |\\
	LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n |\\
	bgzip > "TMP/jeter1.bed.gz"
tabix --comment '#' -f -p bed TMP/jeter1.bed.gz

mv TMP/jeter1.bed.gz "${params.prefix}raw.bed.gz"
mv TMP/jeter1.bed.gz.tbi "${params.prefix}raw.bed.gz.tbi"

cp TMP/jeter2.bed filtered.tsv

sed 's/^chrom/#chrom/' TMP/jeter2.bed |\\
	LC_ALL=C sort -t '\t' -T TMP -k1,1 -k2,2n  |\\
	bgzip > "TMP/jeter2.bed.gz"
tabix --comment '#' -f -p bed TMP/jeter2.bed.gz

mv TMP/jeter2.bed.gz "${params.prefix}filtered.fdr_${FDR}_fold_${fold}.bed.gz"
mv TMP/jeter2.bed.gz.tbi "${params.prefix}filtered.fdr_${FDR}_fold_${fold}.bed.gz.tbi"
"""
}


process PLOT_PEAKS {
tag "${row.chrom}:${row.start}-${row.end} N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(L)
	val(row)
output:
	path("${row.chrom}_${row.start}_${row.end}.pdf"),emit:output
when:
	true
script:
	def chromStart = (row.start as int)
	def chromEnd = (row.end as int)
	def len = 1 + chromEnd - chromStart
	def xStart = java.lang.Math.max(1,chromStart - len)
	def xEnd = java.lang.Math.max(1,chromEnd+len) 
"""
hostname 1>&2
${moduleLoad("ucsc")}
set -o pipefail

mkdir -p TMP

cat << EOF > TMP/sn_bigwig.tsv
${L.join("\n")}
EOF

i=1
cut -f 3 TMP/sn_bigwig.tsv | while read BIGWIG
do
	bigWigToBedGraph -chrom="${row.chrom}" -start=${xStart} -end=${xEnd}  "\${BIGWIG}" stdout |\\
		awk -F '\t' '{X1=int(\$2);X2=int(\$3);for(i=X1;i<X2;++i) {printf("%d\t%s\\n",i+1,\$4);}}' > TMP/depth.\${i}.txt

	i=\$((i+1))
done

${LOAD_DIFFBIND}

cat "${moduleDir}/plot.coverage.R" |\\
	sed 's%__SAMPLESHEET__%TMP/sn_bigwig.tsv%' |\\
	sed 's/__CHROM__/${row.chrom}/g' |\\
	sed 's/__START__/${chromStart}/' |\\
	sed 's/__END__/${chromEnd}/' |\\
	sed 's/__XSTART__/${xStart}/' |\\
	sed 's/__XEND__/${xEnd}/' > TMP/jeter.R

R --vanilla < TMP/jeter.R

mv -v TMP/jeter.pdf "${row.chrom}_${row.start}_${row.end}.pdf"
"""
}


/** merge PDFs using ghostscript */
process MERGE_PDFS {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(L) /* name, list to pdf */
output:
	path("${params.prefix?:""}coverage.pdf"),emit:output
script:
"""
hostname 1>&2
module load gs
mkdir -p TMP

cat << EOF  | awk -F '/' '{printf("%s,%s\\n",\$NF,\$0);}' | LC_ALL=C sort -t, -k1,1V -T TMP | cut -d, -f2  > TMP/jeter.list
${L.join("\n")}
EOF

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=TMP/jeter.pdf @TMP/jeter.list

mv TMP/jeter.pdf "${params.prefix?:""}coverage.pdf"
"""
}



process README {
executor "local"
input:
	path(samplesheet)
	path(occupancy)
	path(affinity)
	path(pca)
	path(ma)
	path(boxplot)
	path(report)
	path(coverage)
output:
	path("${params.prefix}archive.zip"),emit:zip
	path("index.html"),emit:output
script:
	def hg="hg38"
	def nrows=500
"""
cat << __EOF__ > index.html
<html>
<head>
	<title>${params.prefix}</title>
</head>
<body>
<h1>Diffbind output</h1>
<p><a href="https://bioconductor.org/packages/release/bioc/html/DiffBind.html">DiffBind<a>  Compute differentially bound sites
from multiple ChIP-seq experiments using affinity (quantitative) data. Also enables occupancy (overlap) analysis and plotting functions.</p>

<div><quote>Stark R, Brown G (2011). <i>DiffBind: differential binding analysis of ChIP-Seq peak data.</i>
<a href="http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf">http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf</a>.<br/>
(2012). <i>Differential oestrogen receptor binding is associated with clinical outcome in breast cancer.</i>. 
http://www.nature.com/nature/journal/v481/n7381/full/nature10730.html.<br/>
</quote></div>
<h2>Samplesheet</h2>
<table>
<thead>
	<caption>${params.samplesheet}</caption>
</thead>
<tbody>
__EOF__

awk -F '\t' '{tag=(NR==1?"th":"td"); printf("<tr>"); for(i=1;i<=NF;i++) printf("<%s>%s</%s>",tag,\$i,tag);printf("</tr>\\n"); }' ${samplesheet} >> index.html


cat << __EOF__ >> index.html
</tbody>
</table>

<h2>Downsampling</h2>
<p>BAM files were downsampled to the lowest read count per BAM to avoid error : <b>${params.downsample}</b></p>

<h2>Occupancy analysis</h2>
<p>Peaksets provide insight into the potential occupancy of the ChIPed protein at specific genomic regions.</p>
<div><img src="${occupancy.name}"/></div>

<h2>PCA</h2>
<p>To see how well the samples cluster with one another, we can draw a PCA plot</p>
<div><img src="${pca.name}"/></div>

<h2>Affinity analysis</h2>
<p> binding affinity matrix containing a read count for each sample at every consensus binding site, whether or not it was identified as a peak in that sample</p>
<div><img src="${affinity.name}"/></div>

<h2>MA Plot</h2>
<p>MA plots are a useful way to visualize the relationship between the overall binding level at
each site and the magnitude of the change in binding enrichment between conditions, as well
as the effect of normalization on data.</p>
<div><img src="${ma.name}"/></div>

<h2>Box Plot</h2>
<p>Boxplots provide a way to view how read distributions differ between classes of binding sites.</p>
<div><img src="${boxplot.name}"/></div>

<h2>Coverage</h2>
<p>Coverage in WigFiles</p>
<div><a target="pdf"  href="${coverage.name}">${coverage.name}</a></div>

<h2>${nrows} first filtered results.</h2>
<table>
<thead>
        <caption>${params.samplesheet}</caption>
</thead>
<tbody>
__EOF__

gunzip -c *filtered*.bed.gz | head -n ${nrows+1} | awk -F '\t' '{tag=(NR==1?"th":"td"); printf("<tr>"); if(NR==1) {printf("<th>Position</th>");} else {printf("<td><a target=\\"ucsc\\" href=\\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=${hg}&amp;position=%s%%3A%s-%s\\">%s:%s-%s</a></td>",\$1,\$2,\$3,\$1,\$2,\$3);}  for(i=4;i<=NF;i++) printf("<%s>%s</%s>",tag,\$i,tag);printf("</tr>\\n"); }'  >> index.html


cat << __EOF__ >> index.html
</tbody>
</table>
</body>
</html>
__EOF__


mkdir -p "${params.prefix}archive"

cp -v  index.html '${occupancy}' '${affinity}' '${pca}' '${samplesheet}' *.bed.gz *.bed.gz.tbi ${ma} ${boxplot} ${coverage} ${params.prefix}archive/

zip -9r "${params.prefix}archive.zip" "${params.prefix}archive"
"""
}


