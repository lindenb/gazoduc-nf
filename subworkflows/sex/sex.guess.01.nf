/*

Copyright (c) 2023 Pierre Lindenbaum

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

/**
 * Guess sample sex from bam
 *
 */
workflow SEX_GUESS_01 {
	take:
		meta
		sample_bam /* file output of SAMTOOLS_SAMPLES01 */
	main:
		version_ch = Channel.empty()
		sn_ch = sample_bam.splitCsv(header:true,sep:'\t')
		sex_count_ch = SEX_CONTIG_COUNT(meta, sn_ch)
		version_ch = version_ch.mix(sex_count_ch.version)

		sn_sex_ch  = DIGEST(meta, sex_count_ch.output.collect())
		version_ch = version_ch.mix(sn_sex_ch.version)
	emit:
		pdf = sn_sex_ch.pdf
		output = sn_sex_ch.output
		version = version_ch
	}

process SEXUAL_CONTIGS {
executor "local"
input:
	val(meta)
	val(reference)
output:
	path("XY.bed"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2



#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract count from sexual chromosome, guess the sex</entry>
	<entry key="version">${getVersionCmd("awk")}</entry>
</properties>
EOF
"""
}

process SEX_CONTIG_COUNT {
tag "${row.sample}  ${file(row.bam).name}"
cpus 1
afterScript "rm -rf TMP"
input:
	val(meta)
	val(row)
output:
	path("count.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def reference = row.reference?:"NO_FILE"
	def bam = row.bam?:"NO_FILE"
	def mapq = meta.mapq?:30
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

mkdir TMP

test -s "${reference}.fai"
awk -F '\t' '(\$1 ~ /^(chr)?X\$/) {printf("%s\t%s\tX\\n",\$1,\$2);}' "${reference}.fai" > TMP/X.tsv
test -s TMP/X.tsv
awk -F '\t' '(\$1 ~ /^(chr)?Y\$/) {printf("%s\t%s\tY\\n",\$1,\$2);}' "${reference}.fai" > TMP/Y.tsv
test -s TMP/Y.tsv
cat TMP/X.tsv TMP/Y.tsv > TMP/XY.tsv

sort -t '\t' -k1,1 TMP/XY.tsv > TMP/XY.b.tsv
mv TMP/XY.b.tsv TMP/XY.tsv

test `wc -l < TMP/XY.tsv` -eq 2

cut -f 1 TMP/XY.tsv | while read C
do
	echo -ne "\${C}\t" >> TMP/count.txt

	samtools view -@ ${task.cpus} -c  -F 3844 --min-MQ ${mapq} \
		--reference "${reference}" \
		"${bam}" "\${C}" >> TMP/count.txt
done


join -t '\t' -1 1 -2 1 -o '1.1,1.2,1.3,2.2'  TMP/XY.tsv TMP/count.txt |\
	awk -F '\t' '{printf("${row.sample}\t${row.new_sample}\t${bam}\t${reference}\t%s\\n",\$0);}' > count.tsv


#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract count from sexual chromosome, guess the sex</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="bam">${bam}</entry>
	<entry key="reference">${reference}</entry>
	<entry key="version">${getVersionCmd("awk samtools")}</entry>
</properties>
EOF
"""
}

//	awk '{printf("${row.new_sample:?row.sample}\t${bam}\t${reference}\t{contigType}\t{contigLen}\t%s\\n",\$1);}' > count.txt


process DIGEST {
executor "local"
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${meta.prefix?:""}samples.sex.tsv"),emit:output
	path("${meta.prefix?:""}sex.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
	def treshold = meta.treshold?:10
"""
hostname 1>&2
${moduleLoad("R")}
set -o pipefail


cat << EOF > jeter.txt
${L.join("\n")}
EOF


echo -e "sample\tnew_sample\tbam\treference\tfx\tfy\tsex"  > samples.sex.tsv

xargs -a jeter.txt -L 1 cat |\
	sort -T . -t '\t' -k1,1 -k7,7 |\
	awk -F '\t' '(\$7=="X") {FX=int(\$8)/(int(\$6)*1.0);next;} (\$7=="Y") {S="male";FY=int(\$8)/(int(\$6)*1.0);if(FX > (FY * 10.0)) {S="female"};printf("%s\t%s\t%s\t%s\t%f\t%f\t%s\\n",\$1,\$2,\$3,\$4,FX,FY,S);next;}' >> samples.sex.tsv

rm jeter.txt

cat << '__EOF__' | R --vanilla
T1<-read.table("samples.sex.tsv",header = TRUE,sep="\t",comment.char="",stringsAsFactors=FALSE)
male <-T1[T1\$sex=="male",]
head(male)

female <-T1[T1\$sex=="female",]
head(female)
pdf("${meta.prefix?:""}sex.pdf")
plot(1,
	xlab="chrX: count n-read / chrom-length",
	ylab="chrY: count n-read / chrom-length",
	main="Sex guessed from BAMs.",
	sub="${meta.prefix?:""}",
	xlim=c(0,max(T1\$fx)),
	ylim=c(0,max(T1\$fy))
	)

mc <- rgb(0,0,1.0,alpha=0.5)
points(x=male\$fx,y=male\$fy,type='p',col=mc,pch=16)
fc <- rgb(1.0,0,0,alpha=0.5)
points(x=female\$fx,y=female\$fy,type='p',col=fc,pch=16)
legend("topright",legend=c("male","female"),title="Sex",pch=16,col=c(mc,fc)) 
dev.off()
__EOF__

mv 'samples.sex.tsv' '${meta.prefix?:""}samples.sex.tsv'

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract count from sexual chromosome, guess the sex</entry>
	<entry key="treshold">${treshold}</entry>
	<entry key="version">${getVersionCmd("awk R")}</entry>
</properties>
EOF
"""
}
