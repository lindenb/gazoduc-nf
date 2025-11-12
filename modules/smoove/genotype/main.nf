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
process SMOOVE_GENOTYPE {
label "process_short"
tag "${meta.id}"
afterScript  "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(sitesvcf),path(sitesvcfidx)
	tuple val(meta ),path(bam),path(bai)

output:
	tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def sample=meta.id
	def prefix = task.ext.prefix?:meta.id+".genotype"
"""
	hostname 1>&2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP/TMP2
	export TMPDIR=\${PWD}/TMP/TMP2

	smoove genotype -d -x \\
		--vcf "${sitesvcf}" \\
		--outdir TMP \\
		--name ${sample} \\
		--fasta ${fasta} \\
		-p ${task.cpus} \\
		${bam}

	bcftools sort  \\
		--max-mem ${task.memory.giga}G \\
		-T TMP/sort \\
		-o ${prefix}.vcf.gz \\
		-O z \\
		TMP/${sample}-smoove.genotyped.vcf.gz
	
	bcftools index --threads ${task.cpus} -t -f ${prefix}.vcf.gz


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""


stub:
def sample = task.ext.id?:meta.id
"""
touch versions.yml ${sample}-smoove.genotyped.vcf.gz ${sample}-smoove.genotyped.vcf.gz.tbi
"""
}