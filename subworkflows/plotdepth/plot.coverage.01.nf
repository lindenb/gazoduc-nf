/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {DOWNLOAD_GFF3_01} from '../../modules/gff3/download.gff3.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow PLOT_COVERAGE_01 {
	take:
		meta
		reference
		references
		bams
		bed
	main:

		version_ch  = Channel.empty()

		gff_ch = DOWNLOAD_GFF3_01(meta.plus("with_tabix":true),reference)
		version_ch = version_ch.mix(gff_ch.version)


		bams_ch = SAMTOOLS_SAMPLES_01(meta.plus(["with_header":true,"allow_multiple_references":false,"allow_duplicate_samples":true]),reference,references,bams)
		version_ch = version_ch.mix(bams_ch.version)
		
		x_ch = EXTEND_BED(meta,reference,bed)
		version_ch = version_ch.mix(x_ch.version)

		c1_ch = x_ch.bed.splitCsv(header: true,sep:'\t',strip:true)
		c2_ch = bams_ch.output.splitCsv(header: true,sep:'\t',strip:true)

		cov_ch = DRAW_COVERAGE(meta,gff_ch.gff3, c1_ch.combine(c2_ch).map{T->T[0].plus(T[1])})
		version_ch = version_ch.mix(cov_ch.version)

		input_merge_ch = cov_ch.pdf.
			map{T->[T[0].chrom+"_"+T[0].delstart+"_"+T[0].delend, T[0].sample+"\t"+T[1]]}


		merge_pdf_ch = MERGE_PDFS(meta, input_merge_ch.groupTuple())
		version_ch = version_ch.mix(merge_pdf_ch.version)

		version_ch = MERGE_VERSION(meta, "Coverage Plot", "Cov Plotter", version_ch.collect())

		zip_ch = ZIP_ALL(meta,merge_pdf_ch.output.collect())
	emit:
		version = version_ch
		zip = zip_ch.zip
	}


process EXTEND_BED {
	tag "${bed}"
	input:
		val(meta)
		val(reference)
		path(bed)
	output:
		path("extend.bed"),emit:bed
		path("version.xml"),emit:version
	script:
		def extend =meta.extend_bed?:"3.0"
	"""
	hostname 1>&2
	${moduleLoad("jvarkit bedtools")}
	set -o pipefail

	mkdir TMP
	cut -f1,2 "${reference}.fai" > TMP/jeter.genome

	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -R "${reference}" -c 1 "${bed}" | \
		awk -F '\t' '{printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,(NF==3 || \$4==""?".":\$4));}' > TMP/jeter1.bed

	cut -f 1,2,3 TMP/jeter1.bed | bedtools slop -i - -g TMP/jeter.genome ${extend.toString().contains(".")?"-pct":""} -b ${extend} | cut -f 2,3 > TMP/jeter2.bed

	# make title
	echo "chrom\tdelstart\tdelend\ttitle\tstart\tend" > TMP/extend.bed

	# paste to get chrom/original-start/original-end/x- start/x-end
        paste TMP/jeter1.bed TMP/jeter2.bed |\
		sort -T TMP -t '\t' -k1,1V -k2,2n -k3,3n --unique  >> TMP/extend.bed

	test -s TMP/extend.bed

	mv TMP/extend.bed ./

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extends bed</entry>
	<entry key="bedtools.version">\$(bedtools --version)</entry>
</properties>
EOF

	"""
	}


process DRAW_COVERAGE {
	tag "${row.sample} ${row.chrom}:${row.delstart}-${row.delend} ${row.title}"
        afterScript "rm -rf TMP"
	input:
		val(meta)
		path(gtf)
		val(row)
	output:
		tuple val(row),path("${params.prefix?:""}${row.chrom}_${row.start}_${row.end}.${row.sample}.pdf"),emit:pdf
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	${moduleLoad("samtools java R/3.3.3 tabix bedtools")}
	mkdir TMP

	tabix "${gtf.toRealPath()}" "${row.chrom}:${row.start}-${row.end}" |\
		LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		awk -F '\t' '{printf("rect(%d, 0,%d, -1, col=\\\"green\\\")\\n",\$2,\$3);}' > TMP/exons.R

	samtools depth  -a -r "${row.chrom}:${row.start}-${row.end}" "${row.bam}" | cut -f 2,3 > TMP/depth.txt



cat << __EOF__  > TMP/jeter.R

normalize_median <- ${(parseBoolean(params.median))?"TRUE":"FALSE"}


	distance <- 1 + (${row.end} - ${row.start})

	smoothvalue <- 1 + 2 * min((distance-1)%/% 2, ceiling(0.01*distance))

	T<-read.table("TMP/depth.txt", header = FALSE, sep = "\t",colClasses = c("integer","integer"),col.names=c("pos","cov"))
	
	yaxis <- T\\\$cov
	capy <- median(yaxis) *3.0
	
	if( normalize_median == TRUE ) {
		mediany <- median(yaxis)
		if( mediany > 1 ) {
			yaxis <- yaxis / mediany
			}
		}

	## smooth yaxis using runmed
	if(smoothvalue>3) {
		yaxis <- runmed(yaxis,smoothvalue)
		}

	pdf("TMP/jeter.pdf", 15, 5)
	plot(x=T\\\$pos,y=yaxis,
	   main=paste("${row.sample}",ifelse(normalize_median == TRUE ," Normalized","")," Coverage ${row.chrom}:",format(${row.delstart},big.mark=","),"-",format(${row.delend},big.mark=",")," len=",format((${row.delend}-${row.delstart}),big.mark=",")," ","${row.title}",sep=""),
	  sub="${row.bam}",
	  xlab="Position",
          ylab=paste(ifelse(normalize_median == TRUE ," Normalized ",""),"Depth"),
	  xlim = c(min(T[1]),max(T[1])),
	  ylim = c(0,ifelse(normalize_median == TRUE ,2.0,capy)),
	  col=rgb(0.2,0.1,0.5,0.9) , 
	  type="l", 
	  pch=19)
	
	mediany <- ifelse(normalize_median == TRUE ,1.0,median(yaxis))	

	abline(h=mediany, col="blue")
	abline(h=mediany*0.5, col="orange")
	abline(h=mediany*1.5, col="orange")
	abline(v=${row.delstart}, col="green")
	abline(v=${row.delend}, col="green")

	source("TMP/exons.R",local=TRUE)
	dev.off()
__EOF__


R --vanilla < TMP/jeter.R

mv TMP/jeter.pdf "${params.prefix?:""}${row.chrom}_${row.start}_${row.end}.${row.sample}.pdf" 

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">plot coverage</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="bam">${row.bam}</entry>
        <entry key="interval">${row.chrom}:${row.delstart}-${row.delend}</entry>
        <entry key="xinterval">${row.chrom}:${row.start}-${row.end}</entry>
	<entry key="R.version">\$(R --version | head -n1)</entry>
	<entry key="bcftools.version">\$(bcftools version | head -n 2 | paste -s)</entry>
</properties>
EOF
"""
}


process MERGE_PDFS {
	tag "${title} N=${L.size()}"
	input:
		val(meta)
		tuple val(title),val(L)
	output:
		path("${meta.prefix?:""}${title}.pdf"),emit:output
		path("version.xml"),emit:version
	script:
	"""
cat << EOF  | sort -T . -k1,1 -t '\t' | cut -f 2 > jeter.txt
${L.join("\n")}
EOF

	test -s jeter.txt

	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=${meta.prefix?:""}${title}.pdf `cat jeter.txt`

	rm jeter.txt

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extends bed</entry>
        <entry key="gs.version">\$(gs --version)</entry>
        <entry key="count">${L.size()}</entry>
</properties>
EOF
	"""
	}

process ZIP_ALL {
	tag "N=${L.size()}"
	input:
		val(meta)
		val(L)
	output:
		path("${meta.prefix?:""}all.zip"),emit:zip
		path("version.xml"),emit:version
	script:
	"""
hostname 1>&2
	
cat << EOF > jeter.txt
${L.join("\n")}
EOF

zip -9 -j -@ "${meta.prefix?:""}all.zip" < jeter.txt

rm jeter.txt

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">zip pdfs</entry>
</properties>
EOF
	"""
	}
