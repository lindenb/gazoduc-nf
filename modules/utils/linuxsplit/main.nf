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
/** use linux split to split a file into parts */
process LINUX_SPLIT {
label "process_single"
afterScript "rm -rf TMP"
tag "${meta.id?:filein.name}"
input:
        tuple val(meta),path(filein)
output:
        tuple val(meta),path("OUT/*",arity:"0..*"),emit:output
		path("versions.yml"),emit:versions
script:
       	def prefix = task.ext.prefix?:"${meta.id}.chunck."

		if( filein.extension=="gz" && task.ext.suffix==null) {
			throw new IllegalArgumentException("${task.process} when using gz, suffix should be set");
			}
       	def suffix = task.ext.suffix?:".${filein.extension}"
       	def method = task.ext.method?:""
		def filter = task.ext.filter?:""
		def args1 = task.ext.args1?:""
		if(method.trim().isEmpty()) throw new IllegalArgumentException("in ${task.process} method missing");
"""

mkdir -p TMP
${filein.name.endsWith(".gz")?"gunzip -c":"cat"} "${filein}" |\\
	${filter.isEmpty()?"":"${filter}  |"} \\
	split \\
		-a 9 ${args1} \\
		--additional-suffix=${suffix} \\
		${method} \\
		- "TMP/${prefix}"

mv TMP OUT

cat << END_VERSIONS > versions.yml
"${task.process}":
	split: \$(split --version | awk '(NR==1) {print \$NF}')
END_VERSIONS
"""
}
