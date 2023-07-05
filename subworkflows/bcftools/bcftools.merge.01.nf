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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {BCFTOOL_CONCAT_FILE_LIST_01} from '../../modules/bcftools/bcftools.concat.file.list.01.nf'
include {BCFTOOL_CONCAT_COLLECT_01} from '../../modules/bcftools/bcftools.concat.collect.01.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {BCFTOOLS_CONCAT_01} from './bcftools.concat.01.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'

workflow BCFTOOLS_MERGE_01 {
take:
	meta
	bed
	vcfs
main:
	version_ch = Channel.empty()
	d1_ch = SQRT_FILE(meta, vcfs)

	d2_ch = d1_ch.output.splitText().map{it.trim()}
	
	if(bed.name.equals("NO_FILE")) {
		vcf2bed_ch = VCF_TO_BED(meta,vcfs)
		version_ch = version_ch.mix(vcf2bed_ch.version)

		rgns_ch = vcf2bed_ch.chromsbed.splitText(header:false,delim:'\t').
                        map{T->[ T[0], (T[1] as int)+1, T[2] ]}
		}
	else
		{
		rgns_ch = bed.splitText(header:false,delim:'\t').
			map{T->[ T[0], (T[1] as int)+1, T[2] ]}
		}

	d3_ch = MERGE_PART_ONE(meta, d2_ch.combine(rgns_ch))
	version_ch = version_ch.mix(d3_ch.version)

	d4_ch = MERGE_PART_TWO(meta, d3_ch.output.groupTuple())
	version_ch = version_ch.mix(d4_ch.version)

	d5_ch = BCFTOOL_CONCAT_COLLECT_01(meta, d4_ch.vcf.collect() )
	version_ch = version_ch.mix(d5_ch.version)

	 version_ch = MERGE_VERSION(meta, "merge", "merge vcfs", version_ch.collect())
emit:
	vcf = d5_ch.vcf
	index = d5_ch.index
	version = version_ch
}


process MERGE_PART_ONE {
tag "${vcfs} ${contig}:${start}-${end}"
input:
	val(meta)
	tuple val(vcfs),val(contig),val(start),val(end)
output:
	tuple val("${contig}:${start}-${end}"),path("merged.bcf"),emit:output
	path("version.xml")
script:
	def extraBcftoolsMerge = meta.extraBcftoolsMerge?:""
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

bcftools merge ${extraBcftoolsMerge} --regions "${contig}:${start}-${end}" --file-list "${vcfs}" -o merged.bcf -O b

###############################################

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
        <entry key="description">merge variants</entry>
        <entry key="interval">${contig}:${start}-${end}</entry>
        <entry key="vcfs">${vcfs}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}

process MERGE_PART_TWO {
tag "${interval} N=${L.size()}"
afterScript "rm -f jeter.list jeter2.list jeter3.list"
input:
        val(meta)
        tuple val(interval),val(L)
output:
        path("merged.bcf"),emit:vcf
        path("version.xml")
script:
        def extraBcftoolsMerge = meta.extraBcftoolsMerge?:""
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

cat << EOF > jeter.list
${L.join("\n")}
EOF

#
# use the first sample of each VCF to be sure that samples will be ordered the same way
# for the downstream bcftools concat
#
cat jeter.list | while read V
do
        bcftools query -l "\$V" | head -n 1 | tr "\n" "," >> jeter2.list
        echo "\${V}" >> jeter2.list
done

sort -t, -k1,1 jeter2.list | cut -d, -f2 > jeter3.list


bcftools merge ${extraBcftoolsMerge} --regions "${interval}" --file-list "jeter3.list" -o merged.bcf -O b

###############################################

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge variants</entry>
        <entry key="interval">${interval}</entry>
        <entry key="vcfs.count">${L.size()}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}


