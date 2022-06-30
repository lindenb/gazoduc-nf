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
include {moduleLoad;getKeyValue;getModules} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'

BCFTOOLS_CONCAT_01

workflow BCFTOOLS_CONCAT_PER_CONTIG_01 {
	take:
		meta /* params */
		vcfs /* a FILE containing the path to the indexed VCF */
	main:
		version_ch = Channel.empty()

		vcf2bed_ch = VCF_TO_BED(meta,vcfs)
		version_ch = vcf2bed_ch.mix(version_ch.version)

		c1_ch = BED2VCF_PER_CONTIG(meta,vcf2bed_ch.bed)

		c2_ch = c1_ch.output.splitCsv(header:false,sep:"\t")
					
		c3_ch = CONCAT_ONE_CONTIG(meta,c2_ch)

		file_list_ch = COLLECT_TO_FILE_01([:],c3_ch.vcf.collect())
	emit:
		version = vers_ch.version
		
	}

process BED2VCF_PER_CONTIG {
executor "local"
input:
	val(meta)
	path(bed)
output:
	path("vcf2beds.tsv"),emit:output	
script:
"""
hostname 1>&2
${moduleLoad("bedtools")}

cut -f1 "${bed}" | sort | uniq | while read C
do
	awk -vC=\$C -F '\t' '(\$1==C)' "${bed}" | cut -f1,2,3 |\
		sort -T . -t '\t' -k1,1 -k2,2n |\
		bedtools merge > "contig.\${C}.bed"
	awk -vC=\$C -F '\t' '(\$1==C)' "${bed}" | cut -f4 |\
		sort -T . | uniq > "contig.\${C}.vcf.list"

	echo "\$C\t\${PWD}/contig.\${C}.vcf.list\tcontig.\${C}.bed" >> vcf2beds.tsv
done

touch -s vcf2beds.tsv
""
}


process CONCAT_ONE_CONTIG {
input:
	val(meta)
	tuple val(contig),val(vcfs),val(bed)
output:
	path("${meta.prefix?:""}${contig}.merged.bcf"),emit:vcf
	path("${meta.prefix?:""}${contig}.merged.bcf.csi"),emit:index
script:
	def prefix="${meta.prefix?:""}${contig}."
"""
hostname 1>&2
${moduleLoad("bedtools")}
		
	bcftools concat --threads ${task.cpus} \
		--regions-file "${bed}" \
		--no-version --allow-overlaps --remove-duplicates \
		-O b -o "${prefix}merged.bcf" --file-list "${vcfs}"

	bcftools index --threads ${task.cpus}  "${prefix}merged.bcf"
"""
}
