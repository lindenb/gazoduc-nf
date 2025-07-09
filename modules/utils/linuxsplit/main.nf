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
       	def prefix = task.ext.prefix?:"chunck."
       	def suffix = task.ext.suffix?:"."+filein.extension
       	def method = task.ext.method?:""
		if(method.trim().isEmpty()) throw new IllegalArgumentException("in ${task.process} method missing");
"""

mkdir -p TMP
split -a 9 --additional-suffix=${suffix} ${method} "${filein}" "TMP/${prefix}"
mv TMP OUT

cat << END_VERSIONS > versions.yml
"${task.process}":
	split: "todo"
END_VERSIONS
"""
}
