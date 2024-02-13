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

/**
 * Guess sample sex from bam
 *
 */
workflow SEX_GUESS_01 {
	take:
		rows /* contains samples/bam/fasta */
	main:
		version_ch = Channel.empty()
		sex_count_ch = SEX_CONTIG_COUNT(rows)
		version_ch = version_ch.mix(sex_count_ch.version)


		sn_sex_ch  = DIGEST(sex_count_ch.output.collect())

/*
		ch1 = rows.combine(ch0).
			filter{T->T[0].sample.equals(T[1].sample)}.
			map{T->T[0].plus(sex:T[1].sex)}
*/
		version_ch = version_ch.mix(sn_sex_ch.version)
	emit:
		pdf = sn_sex_ch.pdf
//		rows = ch1
		version = version_ch
	}


process SEX_CONTIG_COUNT {
tag "${row.sample}  ${file(row.bam).name}"
cpus 1
afterScript "rm -rf TMP"
input:
	val(row)
output:
	path("samples.sex.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[row.genomeId]
	def reference = genome.fasta
	def bam = row.bam?:"NO_FILE"
	def mapq = params.mapq?:30
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

function countIt {

	awk -F '\t' -vC=\$1 '{if(\$1==C || sprintf("chr%s",C)==\$1) {printf("%s\t0\t%s\\n",\$1,\$2);}}' "${reference}.fai" > "TMP/jeter.bed"
	test -s "TMP/jeter.bed"


	samtools view -@ ${task.cpus} -c  -F 3844 --min-MQ "${mapq}" \\
		-M -L "TMP/jeter.bed" \\
		--reference "${reference}" \\
		"${bam}"  > TMP/count.txt
	test -s TMP/count.txt
	paste TMP/jeter.bed TMP/count.txt | cut -f 1,3,4 > TMP/count2.txt
	}

mkdir -p TMP

test -s "${reference}.fai"

countIt X
mv TMP/count2.txt TMP/countX.txt
countIt Y
mv TMP/count2.txt TMP/countY.txt

paste TMP/countX.txt TMP/countY.txt | awk -F '\t' '{printf("${row.sample}\t${row.bam}\t%s\\n",\$0);}' > samples.sex.tsv



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



process DIGEST {
tag "N=${L.size()}"
input:
	val(L)
output:
	path("${params.prefix}samples.sex.tsv"),emit:output
	path("${params.prefix}sex.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
	def treshold = params.treshold?:10.0
"""
hostname 1>&2
${moduleLoad("R")}
set -o pipefail

echo "sample\tbam\tchrX\tchrX_len\tchrX_count\tchrY\tchrY_len\tchrY_count\tfx\tfy\tsex" > samples.sex.tsv

cat ${L.join(" ")} |\
	sort -T . -t '\t' -k1,1 |\
	awk -F '\t' '{FX=int(\$5)/(int(\$4)*1.0);S="male";FY=int(\$8)/(int(\$7)*1.0);if(FX > (FY * ${treshold} )) {S="female"};printf("%s\t%f\t%f\t%s\\n",\$0,FX,FY,S);next;}' >> samples.sex.tsv



cat << '__EOF__' | R --vanilla
T1<-read.table("samples.sex.tsv",header = TRUE,sep="\t",comment.char="",stringsAsFactors=FALSE)
male <-T1[T1\$sex=="male",]
head(male)

female <-T1[T1\$sex=="female",]
head(female)

pdf("${params.prefix?:""}sex.pdf")
plot(1,
	xlab="chrX: count n-reads / chrom-length",
	ylab="chrY: count n-reads / chrom-length",
	main="Sex guessed from BAMs.",
	sub="${params.prefix?:""}",
	xlim=c(0,max(T1\$fx)),
	ylim=c(0,max(T1\$fy))
	)

mc <- rgb(0,0,1.0,alpha=0.5)
points(x=male\$fx,y=male\$fy,type='p',col=mc,pch=16,
	xlim=c(0,max(T1\$fx)),
	ylim=c(0,max(T1\$fy))
	)
fc <- rgb(1.0,0,0,alpha=0.5)
points(x=female\$fx,y=female\$fy,type='p',col=fc,pch=16,
	xlim=c(0,max(T1\$fx)),
	ylim=c(0,max(T1\$fy))
	)
legend("topright",legend=c("male","female"),title="Sex",pch=16,col=c(mc,fc)) 
dev.off()
__EOF__

mv samples.sex.tsv "${params.prefix}samples.sex.tsv"


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
