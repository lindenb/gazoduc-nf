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

process WALLY_REGION {
tag "${meta.id}"
label "process_single"
label "array100"
conda "${moduleDir}/../../../conda/wally.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(annot_bed),path(annot_tbi)
    tuple val(meta) ,path(bams),path(bai),path(opt_bed)
output:
	tuple val(meta),path("*.png",arity: '1..*'),emit:png
	tuple val(meta),path("*.html",arity: '1..*'),emit:html
	path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def args2 = task.ext.args2?:""
	def att1 = task.ext.att1?:"class='wally'"
	def att2 = task.ext.att2?:"alt=\"${meta.id}\""
	def figcaption = task.ext.figcaption?:"${meta.id}"
	def prefix = task.ext.prefix?:"${meta.id?:"wally"}"
"""
mkdir -p TMP

cat << 'EOF' > TMP/jeter.m4
<figure ${att1}>
   <img ${att2} src="${prefix.isEmpty()?"":"${prefix}."}__PNG__"/>
   <figcaption>${figcaption}</figcatpion>
</figure>
EOF

wally region \\
	${args1} \\
	${args2} \\
	${annot_bed?"-b ${annot_bed}":""} \\
	${opt_bed?"-R \${opt_bed}":""} \\
	-g ${fasta} \\
	${bams}

find . -type f -name "*.png" -printf "%f\\n" | while read P
do
		mv -v "\${P}" "${prefix}.\${P}"
		m4 -P -D__PNG__=\$P < TMP/jeter.m4 > "${prefix}.\${P}.html"

done


cat << EOF > versions.yml
"${task.process}":
	wally: \$(wally --version | grep version -m1)
EOF
"""

stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml "${prefix}.png" 
"""
}

