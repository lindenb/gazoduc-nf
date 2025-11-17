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
process FASTQC {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/fastqc.yml"
input:
	tuple val(meta),path("IN/*")
output:
	tuple val(meta), path("OUT/*.zip", arity: '1..*'),emit:zip /* one zip per fastq */
	tuple val(meta), path("OUT/*.html"),emit:html
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def prefix0 = task.ext.prefix?:""
	def prefix = "${prefix0}${meta.id}."
	def memory_mega = task.memory.mega
	// fastqc doesn't allow more than 10000M
	def fastqc_memory = java.lang.Math.min(memory_mega,10000)
"""
mkdir -p TMP/FQ OUT

find "\${PWD}/IN/" -type l | while read F
do
	ln -v -s "\${F}" "TMP/FQ/${prefix}\${F##*/}"
done

fastqc \\
  --dir TMP \\
  --memory ${fastqc_memory} \\
  --threads ${task.cpus} \\
  -o OUT \\
  -f "fastq" \\
  TMP/FQ/* 1>&2

cat <<-END_VERSIONS > versions.yml
"${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
"""

stub:
"""
touch versions.yml
mkdir -p OUT
touch OUT/${meta.id}.zip
touch OUT/${meta.id}.html
"""
}

