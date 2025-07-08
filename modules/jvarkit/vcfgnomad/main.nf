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

process JVARKIT_VCFGNOMAD {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:vcf.name}"
	input:
		tuple val(meta ),path(gnomadvcf),path(gnomadvcf_idx)
		tuple val(meta ),path(vcf),path(idx)
	output:
		tuple val(meta ),path("*.bcf"),path("*.bcf.csi"),emit:vcf
		path("versions.yml"),emit:versions
	script:
		def buffer_size  = task.ext.buffer_size?:(task.attempt==1?100:10)
		def max_af = task.ext.max_af?:0.01
		def min_af = task.ext.min_af?:0
		def fields = task.ext.fields?:"AF_nfe"
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:""
		def prefix  = task.ext.prefix?:vcf.baseName+".gnomad"
	"""
	set -o pipefail

	mkdir -p TMP

	bcftools view  ${args1} "${vcf}" |\\
	jvarkit -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfgnomad \\
		--bufferSize ${buffer_size} \\
		--min-af ${min_af} \\
		--max-af ${max_af} \\
		--gnomad "${gnomadvcf}" \\
		--fields "${fields}" |\\
	bcftools view ${args2} -O b -o TMP/${prefix}.bcf

	bcftools index --threads ${task.cpus} -f TMP/${prefix}.bcf

	mv TMP/${prefix}.bcf ./
	mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	jvarkit: todo
END_VERSIONS
	"""
	}
