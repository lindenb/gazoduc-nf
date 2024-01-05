/*

Copyright (c) 2024 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'

process SAMTOOLS_DEPTH_PLOT_COVERAGE_01 {
tag "${row.bam}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
	val(reference)
	val(row)
output:
	tuple val(row),path("${meta.prefix?:""}${row.contig}_${row.start}_${row.end}.pdf"),emit:output
	path("version.xml"),emit:version
script:
	def prefix = row.prefix?:""
	def contig = row.contig
	def start = row.start
	def end = row.end
	def mapq = row.mapq?:1

"""
hostname 1>&2
${moduleLoad("bedtools samtools r htslib")}
set -x
mkdir -p TMP

## create BAM list
if ${row.containsKey("bam")} ; then
	echo "${row.bam}" > TMP/bams.list
elif ${row.containsKey("bams")} ; then
	cp "${row.bams}" TMP/bams.list
elif ${row.containsKey("bams_comma")} ; then
	echo "${row.bams_comma}" | tr "," "\\n" | sort | uniq > TMP/bams.list 
else
	echo "no bam source was defined" 1>&2
fi

## reduce BAM list
if ${meta.containsKey("max_bams")} ; then
	head -n '${meta.max_bams}' TMP/bams.list > TMP/bams.list2
	mv TMP/bams.list2 TMP/bams.list
fi

test -s TMP/bams.list

echo "SN\tDEPTH\tmaxDP" > TMP/depths.tsv

# call samtools depth
cat "TMP/bams.list" | samtools samples  | cut -f1,2 | while read SN BAM
do
	
	samtools depth  -a -r "${contig}:${start}-${end}" "\${BAM}" | cut -f 2,3 > TMP/\${SN}.depth.txt
	echo -n "\${SN}\tTMP/\${SN}.depth.txt\t" >> TMP/depths.tsv
	# get max cov
	cut -f2 TMP/\${SN}.depth.txt | sort -T . -n | tail -1 >> TMP/depths.tsv
done

# exons
echo "contig\tstart\tend" > TMP/exons.bed
if ${meta.containsKey("gtf")} ; then
	tabix "${meta.gtf}"  "${contig}:${start}-${end}" |\
	awk -F '\t' '(\$3=="exon")' |\
	cut -f1,4,5 | sort -T TMP | uniq >> TMP/exons.bed
fi


cat << "__EOF__"  > TMP/jeter.R

T2 <- read.table("TMP/depths.tsv", header = TRUE, sep = "\t", colClasses = c("character","character","numeric"),stringsAsFactors=FALSE)
head(T2)
maxDP <- max(T2\$maxDP)
samplesRGB <- rainbow(nrow(T2),alpha=0.3)

# read exon bounds
T3 <- read.table("TMP/exons.bed", header = TRUE, sep = "\t", colClasses = c("character","integer","integer"),stringsAsFactors=FALSE)
head(T3)

pdf("${meta.prefix?:""}${contig}_${start}_${end}.pdf",width=14)

# empty graphics
plot(1, type="l", main="Coverage ${contig}:${start}-${end}",sub="${meta.prefix?:""}",xlab="Genomic position on ${contig}", ylab="DEPTH", xlim=c(${start}, ${end}), ylim=c(0,maxDP))

# plot exons
exonRGB <- rgb(0,0,1.0,alpha=0.05)
for (i in 1:nrow(T3)) {
 rect(xleft=T3[i,]\$start,xright=T3[i,]\$end,ybottom=0,ytop=maxDP,col=exonRGB,border=NA)
 }

# plot each coverage
for (i in 1:nrow(T2)) {
	T1<- read.table(T2[i,]\$DEPTH, header = FALSE, sep = "\t",col.names=c('POS','DP'), colClasses = c("integer","numeric"),stringsAsFactors=FALSE)
	T1\$DP <- runmed(T1\$DP,k=21)
	lines(x=T1\$POS,y=T1\$DP,col=samplesRGB[i])
}

legend("bottomright",legend=T2\$SN,title="Samples",pch=16,col=samplesRGB) 
dev.off()
__EOF__


R --vanilla < TMP/jeter.R


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot coverage with samtools depth</entry>
	<entry key="interval">${contig}:${start}-${end}</entry>
	<entry key="versions">${getVersionCmd("samtools tabix R")}</entry>
</properties>
EOF
"""
}
