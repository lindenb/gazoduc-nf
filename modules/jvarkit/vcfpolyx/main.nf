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

process JVARKIT_VCF_POLYX {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:vcf.name}"
	input:
		tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
		tuple val(meta3),path(dict)
		tuple val(meta ),path(vcf),path(idx)
	output:
		tuple val(meta ),path("*.bcf"),path("*.bcf.csi"),emit:vcf
		path("versions.yml"),emit:versions
	script:
		def size  = task.ext.size?:10
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		def prefix  = task.ext.prefix?:vcf.baseName+".polyx"
		def tag = task.ext.tag?:"POLYX"
	"""
	hostname 1>&2
	set -o pipefail

	mkdir -p TMP

	bcftools view  ${args1} "${vcf}" |\\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfpolyx \\
			-n '${size}' \\
			--reference '${fasta}' \\
			--tag "${tag}"   |\\
	bcftools view ${args2} --write-index -O b -o TMP/${prefix}.bcf

	mv TMP/${prefix}.bcf ./
	mv TMP/${prefix}.bcf.csi ./

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
	"""
	}
