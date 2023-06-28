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
include {getVersionCmd;moduleLoad;parseBoolean} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'

if(!params.containsKey("with_tabix")) throw new IllegalArgumentException("params.with_tabix missing")

/**
	motivation: convert VCF to bed using bedtools query -f '%CHROM\t%POS0\t%END' after split per vcf / split per contig
*/
workflow BCFTOOLS_VCFS_TO_BED_01 {
	take:
		vcf
		bed //restrict to that bed or NO_FILE
	main:
		version_ch = Channel.empty()

		if(vcf.name.endsWith(".list")) {			
			each_vcf = vcf_ch.splitText().map{file(it.trim())}
			rgn_ch = INTERVAL_FOR_VCF(bed,each_vcf)
			}
		else
			{
			rgn_ch = INTERVAL_FOR_VCF(bed,vcf)
			}
		version_ch = version_ch.mix(rgn_ch.version)

		
		each_rgn = rgn_ch.output.splitCsv(header:false,sep:'\t')

		vcfstat_ch = VCF2BED(each_rgn)
		version_ch = version_ch.mix(vcfstat_ch.version)

		beds_ch = MERGE_BEDS(vcfstat_ch.bed.collect())
		version_ch = version_ch.mix(beds_ch.version)

		version_ch = MERGE_VERSION("VCF to BED",version_ch.collect())
	emit:
		version = version_ch
		bed = beds_ch.bed
	}


process INTERVAL_FOR_VCF {
tag "${vcf.name}"
executor "local"
maxForks 5
input:
	path(bed)
	path(vcf)
output:
	path("intervals.tsv"),emit:output
	path("version.xml"),emit:version.xml
script:
"""
hostname 1>&2
${moduleLoad("bcftools bedtools")}
set -o pipefail

bcftools index -s "${vcf.toRealPath()}" |\
	awk -F '\t' '{printf("%s\t0\t%s\t${vcf.toRealPath()}\\n",\$1,\$2);}' |\
	sort -T . -t '\t' -k1,1 -k2,2n > jeter.bed

if ${!bed.name.equals("NO_FILE")} ; then
	bedtools intersect -u -a jeter.bed -b ${bed} |\
		sort -T	. -t '\t' -k1,1	-k2,2n > jeter2.bed
	mv jeter2.bed jeter.bed

fi

	
awk -F '\t' '{printf("%s:%s-%s\t%s\\n",\$1,int(\$2)+1,\$3,\$4);}' jeter.bed > intervals.tsv
rm jeter.bed

###############################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
        <entry key="description">extract intervals from one vcf</entry>
        <entry key="vcf">${vcf.toRealPath()}</entry>
        <entry key="bed">${bed.name}</entry>
        <entry key="versions">${getVersionCmd("bcftools awk bedtools")}</entry>
</properties>
EOF
"""
}

process VCF2BED {
tag "${interval} ${file(vcf).name}"
input:
	val(meta)
	tuple val(interval),val(vcf)
output:
	path("variants.bed"),emit:output
	path("version.xml"),emit:version.xml
script:
"""
hostname 1>&2
${moduleLoad("bcftools bedtools")}
set -o pipefail

bcftools query  -f '%CHROM\t%POS0\t%END' --regions "${interval}"  "${vcf}" |\
	sort -T . -k1,1 -k2,2n -t '\t' |\
	bedtools merge > variants.bed

###############################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
        <entry key="description">extract variants intervals from one vcf</entry>
        <entry key="interval">${interval}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="versions">${getVersionCmd("bcftools bedtools")}</entry>
</properties>
EOF
"""
}

process MERGE_BEDS {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("variants.bed${parseBoolean(meta.with_tabix)?".gz":""}"),emit:output
	path("version.xml"),emit:version.xml
script:
"""
hostname 1>&2
${moduleLoad("bedtools htslib")}
set -o pipefail

sort -T . -k1,1 -k2,2n -t '\t' --merge ${L.join(" ")} |\
	bedtools merge > variants.bed

if ${parseBoolean(meta.with_tabix)} ; then
	bgzip variants.bed
	tabix --force -p bed variants.bed.gz
fi

###############################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
        <entry key="description">merge beds</entry>
        <entry key="bed.count">${L.size()}</entry>
        <entry key="versions">${getVersionCmd("htslib bedtools")}</entry>
</properties>
EOF
"""
}
