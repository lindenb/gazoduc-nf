/*

Copyright (c) 2026 Pierre Lindenbaum

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

process OCTOPUS_TRIOS /* 1 child+2 parents */{
label "process_single"
tag "${meta.id} ${C_bam.name} ${optional_bed?optional_bed.name:""}"
afterScript "rm -rf TMP "
conda "${moduleDir}/../../../conda/octopus.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(C_bam),path(C_bai), // CHILD
	 	path(F_bam),path(F_bai),//PARENT 1
		path(M_bam),path(M_bai), //PARENT 2
		path(optional_bed)
output:
	tuple val(meta),path("*.vcf.gz"),  path("*.vcf.gz.tbi"),path(optional_bed),optional:true,emit:vcf
	path("versions.yml"),emit:versions
script:
	if(!meta.father) throw new IllegalArgumentException("${task.process} missing meta.father")
	if(!meta.mother) throw new IllegalArgumentException("${task.process} missing meta.mother")
	if(meta.father.equals(meta.mother)) throw new IllegalArgumentException("${task.process} father=mother")
	if(meta.father.equals(meta.id)) throw new IllegalArgumentException("${task.process} father=id")
	if(meta.mother.equals(meta.id)) throw new IllegalArgumentException("${task.process} mother=id")
	def child = meta.id
	def father = meta.father
	def mother = meta.mother
	def error_model = task.ext.error_model?:"PCRF.NOVASEQ"
	def prefix = task.ext.prefix?:child
	def args1 = task.ext.args1?:""
"""
	hostname 1>&2
	mkdir -p TMP/TMP 
	export TMPDIR=\${PWD}/TMP/TMP	


	## C'EST MORT https://github.com/luntergroup/octopus/issues/39

	octopus \\
		${args1} \\
		--working-directory TMP/TMP \\
		-R "${fasta}" \\
		--reads "${C_bam}"  "${F_bam}"  "${M_bam}" \\
		--sequence-error-model ${error_model} \\
		--paternal-sample "${father}" \\
	    --maternal-sample "${mother}" \\
		--threads ${task.cpus} \\
		${optional_bed?"--regions-file ${optional_bed}":""} \\
		-o TMP/jeter.vcf.gz \\
		1>&2


	find TMP -type f 1>&2

	mv TMP/jeter.vcf.gz        "${prefix}.vcf.gz"
	mv TMP/jeter.vcf.gz.tbi    "${prefix}.vcf.gz.tbi"


cat << EOF > versions.yml
"${task.process}":
    octopus: todo
	error_model: ${error_model}
EOF
"""
}
