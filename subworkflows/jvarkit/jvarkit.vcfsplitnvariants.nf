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
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf' addParams(with_header:false)
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'


if(!params.containsKey("split_vcf_method")) throw new IllegalArgumentException("split_vcf_method missing")
if(!params.containsKey("extraBcfTools")) throw new IllegalArgumentException("extraBcfTools missing")

/**
 * call jvarkit vcfsplitnvariants one or more vcf
 *
 */
workflow JVARKIT_VCF_SPLIT_N_VARIANTS_01 {
	take:
		vcf /* path to vcf, or file with .list suffix */
		bed /* limit to that BED or NO_FILE */
	main:
		version_ch = Channel.empty()
		bed_ch = VCF_TO_BED(vcf)
		version_ch = version_ch.mix(bed_ch.version)

		contig_vcf_ch = bed_ch.bed.splitCsv(sep:'\t',header:false).
			map{T->[T[0],T[3]]}

		rch = SPLIT_VCF(bed, contig_vcf_ch)
		version_ch = version_ch.mix(rch.version)

		concat_ch = CONCAT_FILES_01([:],rch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		version_ch = MERGE_VERSION("jvarkit2intervals",version_ch.collect())

	emit:
		output = concat_ch.output
		version = version_ch
	}



process SPLIT_VCF {
tag "${contig} ${file(vcf).name}"
afterScript "rm -rf TMP"
input:
	path(userbed)
	tuple val(contig),val(vcf)
output:
	path("vcfs.list"),emit:output
	path("version.xml"),emit:version
script:
	def method = params.split_vcf_method ?:" --vcf-count  1000"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit bedtools")}
set -o pipefail
mkdir -p OUT TMP

if ${!userbed.name.equals("NO_FILE")} ; then
	${userbed.name.endsWith(".gz")?"gunzip -c":"cat"} "${userbed}" |\
		awk -F '\t' '(\$1=="${contig}")' |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge > TMP/jeter.bed

	# prevent empty bed
	if test ! -s TMP/jeter.bed ; then
		echo -e 'K_\t0\t1' > TMP/jeter.bed
	fi
fi

bcftools view ${params.extraBcfTools?:""} ${userbed.name.equals("NO_FILE")?"":" --regions-file TMP/jeter.bed "} "${vcf}"  ${userbed.name.equals("NO_FILE")?"\"${contig}\"":""} |\
	java -jar \${JVARKIT_DIST}/vcfsplitnvariants.jar ${method} -o \${PWD}/OUT/split.${contig}

find OUT/ -type f -name "*.vcf.gz" | while read F
do
	bcftools index -t --force "\${F}"
done

find \${PWD}/OUT/ -type f -name "*.vcf.gz" > vcfs.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split VCF per contig and using jvarkit/vcfsplitnvariants</entry>
        <entry key="method">${method}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="bed">${userbed}</entry>
        <entry key="contig">${contig}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/vcfsplitnvariants")}</entry>
</properties>
EOF
"""
}

