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

process REGENIE_MAKE_ANNOT {
tag "${meta1.id} ${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"   
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(annotations)
	tuple val(meta ),path(tsv)
output:
    tuple val(meta),path("OUT/*.aaf.txt.gz", arity: '0..*'),optional:true,emit:aaf
    tuple val(meta),path("OUT/*.annot.txt.gz", arity: '0..*'),optional:true,emit:annot
    tuple val(meta),path("OUT/*.mask.txt.gz", arity: '0..*'),optional:true,emit:mask
    tuple val(meta),path("OUT/*.setfile.txt.gz", arity: '0..*'),optional:true,emit:setfile
	path("versions.yml"),emit:versions
script:
	def jvm = task.ext.jvm?:"-Djava.io.tmpdir=TMP "
	def jvarkit = task.ext.jvarkit?:"java -jar  ${jvm} \${HOME}/jvarkit.jar"
    def prefix = task.ext.prefix?:"${meta.id}"
	def args1 = task.ext.args1?:""
"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT/

${tsv.name.endsWith(".gz")?"gunzip -c":"cat"} "${tsv}" |\\
${jvarkit} regeniemakeannot \\
		-m "${annotations}" \\
		--prefix "${prefix}." \\
		--reserve 20 \\
		-o \${PWD}/OUT \\
		${args1} \\
		--gzip \\
		-N 5000

find OUT -type f 1>&2


cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(${jvarkit} --version)"
EOF
"""
stub:
"""
 def prefix = task.ext.prefix?:"${meta.id}"
touch versions.yml ${prefix}.tsv.gz
"""
}