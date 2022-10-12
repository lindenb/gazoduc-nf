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

include { VCF_TO_BED } from '../../modules/bcftools/vcf2bed.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {md5;getVersionCmd;moduleLoad;isHg38;isHg19} from '../../modules/utils/functions.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'

workflow VCF_RETROCOPY_01 {
	take:
		meta	
		reference
		vcf
		gtf
		gff3
		bams
	main:
		version_ch = Channel.empty()
		ch1_ch  = VCF_TO_BED(meta,vcf)
		version_ch = version_ch.mix(ch1_ch.version)

		ch2_ch = ch1_ch.bed.splitCsv(header:false,sep:'\t').
			map{T->[T[0],T[3]]}

		ch3_ch = DOWNLOAD_KNOWN(meta,reference,gtf)
		version_ch = version_ch.mix(ch3_ch.version)

		ch4_ch = SCAN_RETROCOPY(meta,reference,ch3_ch.output,gtf,ch2_ch)
		version_ch = version_ch.mix(ch4_ch.version)
		


		ch5_ch = COLLECT_TO_FILE_01(meta, ch4_ch.vcf.collect())
		version_ch = version_ch.mix(ch5_ch.version)

		ch6_ch = BCFTOOLS_CONCAT_01(meta,ch5_ch.output)
		version_ch = version_ch.mix(ch6_ch.version)


		if(!bams.name.equals("NO_FILE")) {
			samples_ch = SAMTOOLS_SAMPLES_01(
				meta.plus("with_header":false,"allow_multiple_references":false),
				reference,
				file("NO_FILE"),
				bams)
			version_ch = version_ch.mix(samples_ch.version)

			join_ch=JOIN_SAMPLE_BAM(meta,reference,ch6_ch.vcf, samples_ch.output)
			version_ch = version_ch.mix(join_ch.version)

			plot_ch = PLOT(meta,reference,gff3, join_ch.bed.splitCsv(header:true,sep:'\t'))
			version_ch = version_ch.mix(plot_ch.version)

			zip_ch = ZIP_PLOT(meta,plot_ch.output.collect()).zip
			}
		else
			{
			zip_ch = file("NO_FILE")
			}

		version_ch = MERGE_VERSION(meta, "VCF retrocopies", "VCF retrocopies", version_ch.collect())
	emit:
		vcf = ch6_ch.vcf
		version = version_ch
		zip = zip_ch
	}

// cf. https://gist.github.com/lindenb/3ee8628039713f5baa532e6f0d6454df
process DOWNLOAD_KNOWN {
input:
	val(meta)
	val(reference)
	val(gtf)
output:
	path("known.txt"),emit:output
	path("names.txt"),emit:gene_names
	path("version.xml"),emit:version
script:
	def agent="Mozilla/5.0 (X11; Linux i686; rv:103.0) Gecko/20100101 Firefox/103.0"
	def url ="http://retrogenedb.amu.edu.pl/static/download/homo_sapiens.tsv"
	def base="https://www.bioinfo.mochsl.org.br"
if(isHg38(reference) || isHg19(reference))
"""
hostname 1>&2
${moduleLoad("jvarkit")}
set -o pipefail

java -jar \${JVARKIT_DIST}/gtf2bed.jar "${gtf}" -c gtf.feature,gene_name,gene_id | cut -f 4,5,6 |\
awk -F '\t' '\$1=="gene"' | cut -f 2,3 |\
sort -t '\t' -T . -k1,1 | uniq > jeter.gene.gene_id.txt


wget -O - "${url}" | cut -f 3 |\
	tail -n +2 | awk '(\$1!="." && \$1!="-" && \$1!="")' > tmp1.txt

wget  -U "${agent}" -O jeter.html \
	--keep-session-cookie --save-cookies=cookies.txt \
	"${base}/rcpedia/"
test -s cookies.txt


wget -O jeter.html  \
	-U "${agent}" \
	--keep-session-cookie --save-cookies=cookies.txt \
	--load-cookies=cookies.txt \
	"${base}/rcpedia/species/view/1"

wget -O jeter.html  \
	-U "${agent}" \
	--load-cookies=cookies.txt \
	--keep-session-cookie --save-cookies=cookies.txt \
	--post-data="data%5BRetrocopies%5D%5Bspecie_id%5D=1&data%5BRetrocopies%5D%5Bsearch_string%5D=" \
	"${base}/rcpedia/retrocopies/search"

##  number of pages should be checked in the future
for P in {1..157}
do 
	wget \
		-U "${agent}" \
		--keep-session-cookie  --save-cookies=cookies.txt --load-cookies=cookies.txt \
		-O jeter.html "${base}/rcpedia/retrocopies/search/page:\${P}"
	
	cat jeter.html | tr '"' '\n' | grep searchbygenename |\
		awk -F '/' '{print \$NF}' | sort | uniq > tmp2.txt
	paste -s -d, < tmp2.txt
	cat tmp2.txt >> tmp1.txt
	rm tmp2.txt
done

cat tmp1.txt | sort | uniq >  names.txt
rm tmp1.txt

join -t '\t' -1 1 -2 1 -o '2.2' names.txt jeter.gene.gene_id.txt |\
sort | uniq > known.txt

rm jeter.gene.gene_id.txt

##### 
cat << EOF > version.xml
<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">download known retrogenes</entry>
		<entry key="base"><a>${base}</a></entry>
		<entry key="url"><a>${url}</a></entry>
		<entry key="count">\$( wc -l known.txt)</entry>
</properties>
EOF
"""
else
"""
cat << EOF > version.xml
<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">download known retrogenes: undefined. not human</entry>
		<entry key="count">0</entry>
</properties>
EOF
touch "known.txt"
"""
}

process SCAN_RETROCOPY {
tag "${contig} ${vcf}"
afterScript "rm -rf TMP"
memory "5g"
input:
	val(meta)
	val(reference)
	val(known)
	val(gtf)
	tuple val(contig),val(vcf)
output:
	path("${contig}.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail
${moduleLoad("bcftools jvarkit htslib")}

mkdir TMP
tabix "${gtf}" "${contig}" > TMP/jeter.gtf

bcftools view -O z -o TMP/jeter.vcf.gz "${vcf}" "${contig}"
bcftools index --tbi TMP/jeter.vcf.gz

java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  -jar ${JVARKIT_DIST}/gtfretrocopy.jar \
		--gtf "TMP/jeter.gtf" \
		--known "${known}" \
		TMP/jeter.vcf.gz > TMP/jeter.vcf

bcftools sort --max-mem "${task.memory.giga}G" -T tmp -o "${contig}.bcf" -O b TMP/jeter.vcf
bcftools index "${contig}.bcf"

###
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">scan vcf for retrogenes</entry>
	<entry key="contig">${contig}</entry>
	<entry key="versions">${getVersionCmd("jvarkit/gtfretrocopy bcftools tabix")}</entry>
</properties>
EOF

"""
}

process JOIN_SAMPLE_BAM {
executor "local"
input:
	val(meta)
	val(reference)
	path(vcf)
	path(samples)
output:
	path("output.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools bedtools")}
set -o pipefail

cut -f 1,3 "${samples}" |\
	sort -T . -t '\t' -k1,1 --unique > jeter1.txt

bcftools query -f '[%SAMPLE\t%CHROM\t%POS\t%END\t%GT\\n]' "${vcf}" |\
	grep -v '0/0\$' |\
	grep -v '\\./\\.\$' |\
	sort -T . -t '\t' -k1,1  > jeter2.txt

echo -e "contig\tstart\tend\tbams" > output.bed

join -t '\t' -1 1 -2 1 -o '2.2,2.3,2.4,1.2' jeter1.txt jeter2.txt |\
	sort -t '\t' -T . -k1,1 -k2,2n |\
	bedtools merge -c 4 -o distinct >> output.bed


###
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">join samples and VCF</entry>
	<entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""	
}

process PLOT {
tag "${row.contig}:${row.start}-${row.end}"
memory "5g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(gff3)
	val(row)
output:
	path("${meta.prefix?:""}${row.contig}_${row.start}_${row.end}.html"),emit:output
	path("version.xml"),emit:version
script:
	def mapq = meta.mapq?:"30"
	def uid = md5("${row.contig}:${row.start}-${row.end}")
"""
hostname 1>&2
${moduleLoad("jvarkit")}

mkdir -p TMP

echo '${row.bams}' | tr "," "\\n" | head -n 10 > TMP/jeter.list

java  -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP  -jar ${JVARKIT_DIST}/coverageplotter.jar \
	-R "${reference}" \
	--mapq "${mapq}" \
	--include-center \
	${gff3.name.equals("NO_FILE")?"":"--gff3 '${gff3.toRealPath()}'"} \
	--region "${row.contig}:${row.start}-${row.end}" \
	TMP/jeter.list > TMP/jeter.html


cat << EOF > TMP/jeter.xsl
<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:u="https://umr1087.univ-nantes.fr/" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:svg="http://www.w3.org/2000/svg" version="1.0">
  <xsl:output method="xml" omit-xml-declaration = "yes"/>
  <xsl:template match="svg:metadata">
  <svg:metadata>
    <xsl:apply-templates select="@*|node()" /> 
            <u:Variant rdf:about="${uid}">
              <u:build>hs37d5_all_chr</u:build>
              <u:contig>${row.contig}</u:contig>
              <u:start>${row.start}</u:start>
              <u:end>${row.end}</u:end>
              <u:id>${uid}</u:id>
              <u:filter>PASS</u:filter>
              <u:svtype>RETROCOPY</u:svtype>
              <u:svlen>${1 + (row.end as int) - (row.start as int)}</u:svlen>
            </u:Variant>
  </svg:metadata>
  </xsl:template>

  <xsl:template match="div[@id='__PLACEHOLDER__']">
   <div>
   <h3>Description</h3>
   <div></div>
   </div>
  </xsl:template>

  <xsl:template match="@*|node()">
    <xsl:copy>
      <xsl:apply-templates select="@*|node()" />
    </xsl:copy>
  </xsl:template>

</xsl:stylesheet>
EOF


xsltproc TMP/jeter.xsl TMP/jeter.html >  "${meta.prefix?:""}${row.contig}_${row.start}_${row.end}.html"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot Retrogene</entry>
	<entry key="interval">${row.contig}:${row.start}-${row.end}</entry>
	<entry key="mapq">${mapq}</entry>
	<entry key="gff3">${gff3}</entry>
	<entry key="versions">${getVersionCmd("jvarkit/coverageplotter xsltproc")}</entry>
</properties>
EOF
"""
}

process ZIP_PLOT {
executor "local"
input:
	val(meta)
	val(L)
output:
	path("${meta.prefix?:""}plots.zip"),emit:zip
script:
"""
cat << EOF > jeter.list
${L.join("\n")}
EOF

zip -j -@ "${meta.prefix?:""}plots.zip" < jeter.list
"""
}
