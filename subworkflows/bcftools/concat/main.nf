include {BCFTOOLS_CONCAT_CONTIGS} from '../concat.contigs'

def ALL_CONTIGS="ALL"


workflow BCFTOOLS_CONCAT {
take:
	vcfs // tuple [vcf,idx]
	bed //optional bed file
main:
	ch1 = CONTIGS(vcfs,bed)
	ch2 = ch1.output.splitText().
		map{[it[0].trim(),it[1],it[2]]}

	ch4 = BCFTOOLS_CONCAT_CONTIGS(ch2, bed)
emit:	
	output = ch4.output

}

process CONTIGS {
tag "${vcf.name}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
label "process_single"
afterScript "rm -rf TMP"
input:
	tuple path(vcf),path(idx)
	path(bed)
output:
	tuple path("chroms.txt"),path(vcf),path(idx),emit:output
script:
	def method = task.ext.method?:"per_contig"

if(method.equals("per_contig") && !bed.name.contains("."))
"""
set -o pipefail
bcftools index -s "${vcf}" | cut -f1 > "chroms.txt"
"""
else if(method.equals("per_contig"))
"""
set -o pipefail
cut -f1 '${bed}' | uniq | sort | uniq > a.txt
bcftools index -s "${vcf}" | cut -f1 | sort | uniq > b.txt
comm -12 a.txt b.txt > "chroms.txt"
rm a.txt b.txt
"""
else if(method.equals("all"))
"""
echo "${ALL_CONTIGS}" > chroms.txt
"""
else
"""
echo "not a valid method ${method}"
"""
}
