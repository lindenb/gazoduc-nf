/*

Copyright (c) 2026 Pierre Lindenbaum

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

/**
 * Guess sample sex from bam
 *
 */
workflow SEX_GUESS {
	take:
		meta
		fasta
		fai // channel ! [meta,fai]
		dict
		bam // [[id:id],bam,bai]
	main:
		versions = Channel.empty()

		chrXY_ch = fai.splitCsv(sep:'\t',header:false)
			.filter{it[0].matches("chr[XY]")}
			.map{[it[0],it[1]]}//chrom,length


		sex_count_ch = SEX_CONTIG_COUNT(
			fasta,
			fai,
			bam.combine(chrXY_ch)
			)
		versions = versions.mix(SEX_CONTIG_COUNT.out.versions)

		SEX_CONTIG_COUNT.out.tsv
			.map{it[1]}
			.splitCsv(sep:'\t',header:false)
			branch{v->
				KX: v[1].matches("(chr)?X")
				KY: v[1].matches("(chr)?Y")
			}

		sn_sex_ch  = DIGEST(sex_count_ch.output.collect())

/*
		ch1 = rows.combine(ch0).
			filter{T->T[0].sample.equals(T[1].sample)}.
			map{T->T[0].plus(sex:T[1].sex)}
*/
		version_ch = version_ch.mix(sn_sex_ch.version)
	emit:
		pdf = sn_sex_ch.pdf
		versions = versions
	}


process SEX_CONTIG_COUNT {
label "process_single"
tag "${meta.id?:} ${contig}"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path(bam),path(bai),val(contig),val(chromLen)
output:
	tuple val(contig),path("*.tsv"),emit:tsv
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def prefix = task.ext.prefix?:"${meta.id}.${contig}"
"""
# get coverage on this contig
samtools coverage ${args1} \\
	--reference "${fasta}" \\
	-r "${contig}:1-${chromLen}" \\
	--no-header "${bam}" |\\
	awk -F '\t' '{printf("${meta.id}\t${contig}\t${chromLen}\t%s\\n",\$7);}' > TMP/depth.txt

mv TMP/depth.txt ${prefix}.tsv

cat << EOF > versions.yml
${task.process}:
	samtools: todo
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

echo "sample\tbam\tchrX\tdpx\tchrY\tdpy\tsex" > samples.sex.tsv

cat ${L.join(" ")} |\
	sort -T . -t '\t' -k1,1 |\
	awk -F '\t' '{FX=\$4;FY=\$6;if(FX > (FY * ${treshold} )) {S="female"};printf("%s\t%s\\n",\$0,S);next;}' >> samples.sex.tsv



cat << '__EOF__' | R --vanilla
T1<-read.table("samples.sex.tsv",header = TRUE,sep="\t",comment.char="",stringsAsFactors=FALSE)
male <-T1[T1\$sex=="male",]
head(male)

female <-T1[T1\$sex=="female",]
head(female)

pdf("${params.prefix?:""}sex.pdf")
plot(1,
	xlab="Depth chrX",
	ylab="Depth chrY",
	main="Sex guessed from BAMs.",
	sub="${params.prefix?:""}",
	xlim=c(0,max(T1\$dpx)),
	ylim=c(0,max(T1\$dpy))
	)

mc <- rgb(0,0,1.0,alpha=0.5)
points(x=male\$dpx,y=male\$dpy,type='p',col=mc,pch=16,
	xlim=c(0,max(T1\$dpx)),
	ylim=c(0,max(T1\$dpy))
	)
fc <- rgb(1.0,0,0,alpha=0.5)
points(x=female\$dpx,y=female\$dpy,type='p',col=fc,pch=16,
	xlim=c(0,max(T1\$dpx)),
	ylim=c(0,max(T1\$dpy))
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
