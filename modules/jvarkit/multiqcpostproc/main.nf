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

process JVARKIT_MULTIQCPOSTPROC {
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	tag "${meta.id?:vcf.name}"
	input:
		tuple val(meta1),path(sample2pop)
		tuple val(meta ),path(multiqc_data_dir)
	output:
		tuple val(meta ),path("*.json",arity:"0..*"),emit:json
		path("versions.yml"),emit:versions
	script:
		def jvm= task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
mkdir -P TMP/OUT
jvarkit ${jvm} \\
	multiqcpostproc \\
	--sample2collection "${sample2pop}" \\
	-o TMP/OUT \\
	${multiqc_data_dir}/

if ls TMP/OUT/*.json
then
	mv TMP/OUT/*.json ./
fi

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""
	
stub:
"""
touch versions.yml
touch ${meta.id}.polyx.bcf
touch ${meta.id}.polyx.bcf.csi
"""		
}
