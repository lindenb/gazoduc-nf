/* author Pierre Lindenbaum */

params.vcfs="NO_FILE"
params.samples="NO_FILE"

workflow {
	each_vcf = Channel.fromPath(params.vcfs).splitText().map{it.trim()}

	c2vcf = CHROMS_IN_VCF(each_vcf)

	xch = EXTRACT_SAMPLE(file(params.samples), c2vcf.out.splitCsv(header:false,sep:',') )

	CONCAT(xch.output.groupTuple())
	}

process CHROMS_IN_VCF {
executor "local"
input:
	val(vcf)
output:
	path("chroms.tsv"),emit:output
script:
"""
set -o pipefail
bcftools stats -s "${vcf}" | awk -F '\t' '{prinft("%s,${vcf}\\n",\$1);}' > chroms.tsv
"""
}


process EXTRACT_SAMPLE {
tag "${contig} / ${vcf}"
cpus 6
input:
	path(samples)
	tuple val(contig),val(vcf)
output:
	tuple val(vcf),path("contig.bcf"),emit:output
script:
"""
set -o pipefail
bcftools view --threads 5 -O u -S "${samples}" "${vcf}" "${contig}" | bcftools view --min-ac 1 -O -o b contig.bcf
bcftools index --threads ${task.cpus} contig.bcf
"""
}


process CONCAT {
tag "${vcf} N=${L.size()}"
cpus 6
input:
	tuple val(vcf),val(L)
output:
	tuple path("${file(vcf).getSimpleName()}.bcf"),emit:output
script:
"""
cat << EOF > tmp.list
${L.join("\n")}
EOF
bcftools concat --threads ${task.cpus} --file-list tmp.list -a -o "${file(vcf).getSimpleName()}.bcf"
bcftools index --threads ${task.cpus} "${file(vcf).getSimpleName()}.bcf"
rm tmp.list
"""
}
