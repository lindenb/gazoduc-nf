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
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from './jvarkit.vcf2intervals.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'


/**
 * call call JVARKIT_VCF_TO_INTERVALS_01 and group by cluster
 *
 */
workflow JVARKIT_VCF_TO_BED_01 {
	take:
		meta
		reference /* not used, just in case.. */
		vcf /* path to vcf, or file with .list suffix */
		userbed /* limit to that BED or NO_FILE (exon.bed) */
	main:
		version_ch = Channel.empty()
		if(userbed.name.equals("NO_FILE")) {
			ch1_ch = JVARKIT_VCF_TO_INTERVALS_01(meta, reference, vcf, file("NO_FILE"))
			version_ch = version_ch.mix(ch1_ch.version)

			interval_vcf_ch = ch1_ch.bed	
			}
		else
			{
			bed_ch = VCF_TO_BED(meta, vcf)
			version_ch = version_ch.mix(bed_ch.version)

			ch2_ch = INTERSECT_BED(meta, bed_ch.bed, userbed)
			version_ch = version_ch.mix(ch2_ch.version)			

			interval_vcf_ch = ch2_ch.bed
			}

		each_contig_vcf_ch = interval_vcf_ch.splitCsv(header:false,sep:'\t').map{T->[T[0],T[3]]}.unique()

		ch2_ch = CLUSTER_BED(meta, reference, interval_vcf_ch , each_contig_vcf_ch)
		version_ch = version_ch.mix(ch2_ch.version)
		
		concat_ch = CONCAT_FILES_01(meta, ch2_ch.output.collect())

		version_ch = MERGE_VERSION(meta,"jvarkit2bed","VCF to bed",version_ch.collect())
	emit:
		output = concat_ch.output /* path to bed , path to vcf */
		version = version_ch
	}


process INTERSECT_BED {
tag "${vcfbed.name} / ${userbed.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(vcfbed)
	path(userbed)
output:
	path("intersect.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools")}
set -o pipefail

mkdir -p TMP

sort -t '\t' -T TMP -k1,1 -k2,2n "${vcfbed}" | cut -f1-4 >  TMP/a.bed
sort -t '\t' -T TMP -k1,1 -k2,2n "${userbed}" | cut -f1,2,3 | bedtools merge >  TMP/b.bed

# 4th column contains the path to the VCF
bedtools intersect -a TMP/a.bed -b TMP/b.bed |\
	sort -t '\t' -T TMP -k1,1 -k2,2n > intersect.bed

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">intersection between vcf bed and user bed</entry>
	<entry key="version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}



process CLUSTER_BED {
tag "${contig}/${bed}/${vcf}"
memory "5g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
	tuple val(contig),val(vcf)
output:
	path("beds.${contig}.list"),emit:output
	path("version.xml"),emit:version
script:
	def method = meta.jvarkit_bedcluster_method?:"--size '1mb'"
"""
hostname 1>&2
${moduleLoad("bedtools jvarkit")}
set -o pipefail

mkdir TMP BEDS

test ! -z "${method}"

awk -F '\t' '\$1=="${contig}" && \$4=="${vcf}"' "${bed}" |\
	cut -f 1,2,3 |\
	sort -t '\t' -T TMP -k1,1 -k2,2n |\
	bedtools merge |\
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedcluster.jar \
		-R "${reference}" \
		${method} \
		-o BEDS 

find \${PWD}/BEDS -type f -name "*.bed" | awk '{printf("%s\t${vcf}\\n",\$0);}' > beds.${contig}.list
test -s beds.${contig}.list


#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">intersection between vcf bed and user bed</entry>
	<entry key="contig">${contig}</entry>
	<entry key="method">${method}</entry>	
	<entry key="vcf">${vcf}</entry>	
	<entry key="version">${getVersionCmd("bedtools awk jvarkit/bedcluster")}</entry>
</properties>
EOF
"""
}


