/*

Copyright (c) 2025 Pierre Lindenbaum

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
process XHUNTER_APPLY {
	label "process_single"
	tag "${meta.id}"
	afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/xhunter.yml"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(catalog)
		tuple val(meta ),path(bam),path(bai)
	output:
		tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
		tuple val(meta),path("*.json"),emit:json
		tuple val(meta),path("*.bam"),path("*.bai"),optional:true,emit:bam
		path("versions.yml"),emit:versions
	script:
		def keepbam = ((task.ext.keep_bam?:false) as boolean)
		def mode  = (task.ext.mode?:"streaming")
		def prefix = task.ext.prefix?:"${meta.id}"
		def sex= meta.sex?:"undefined"
		def args1 = task.ext.args1?:"--region-extension-length  1000 "
	"""
	hostname 1>&2
	mkdir -p TMP

	ExpansionHunter \\
		${args1} \\
		--reads "${bam}" \\
		--reference "${fasta}" \\
		--variant-catalog "${catalog}" \\
		--threads ${task.cpus} \\
		--output-prefix "TMP/jeter" \\
		--analysis-mode "${mode}" \\
		--sex '${sex}' \\
		--log-level info 1>&2
	
	# rename sample
	bcftools query -l TMP/jeter.vcf | awk '{printf("%s\t${meta.id}\\n",\$1);}' > TMP/rename.txt

	bcftools reheader --fai "${fai}" --samples TMP/rename.txt TMP/jeter.vcf |\\
		awk '/^#CHROM/ {printf("##samples.sex=${meta.id}=${meta.sex}\\n",S);} {print;}' > TMP/jeter2.vcf
	
	mv TMP/jeter2.vcf TMP/jeter.vcf
	
	#sort and index
	bcftools annotate --set-id '%INFO/REPID' -O u TMP/jeter.vcf |\\
		bcftools sort --max-mem '${task.memory.giga}G' -T TMP/sort -O z -o "TMP/${prefix}.vcf.gz"
	bcftools index   --threads ${task.cpus} -f -t "TMP/${prefix}.vcf.gz"


	if  ${keepbam} ; then
		samtools sort -T TMP/sort  --threads ${task.cpus} -O BAM --reference "${fasta}" -o "${prefix}.realigned.cram" "TMP/jeter_realigned.bam"
		samtools index --threads ${task.cpus} "${prefix}.realigned.bam"
	fi


	mv TMP/${prefix}.vcf.gz ./
	mv TMP/${prefix}.vcf.gz.tbi ./
	mv TMP/jeter.json "${prefix}.json"

cat << EOF > versions.yml
"${task.process}":
    xhunter: \$( ExpansionHunter --version 2>&1 | awk '/^Start/ {print \$NF}')
EOF
"""
stub:
"""
touch versions.yml "${meta.id}.json" "${meta.id}.vcf.gz" "${meta.id}.vcf.gz.tbi" "${meta.id}.realigned.bam" "${meta.id}.realigned.bam.bai"
"""
}
