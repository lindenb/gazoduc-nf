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
process SEQKIT_SPLIT2 {
label "process_single"
conda "${moduleDir}/../../../conda/seqkit.yml"
tag "${meta.id} ${R1.name} ${R2?R2.name:""}"
afterScript "rm -rf TMP"
// cpus 4 set in config file
input:
	val(meta),path(R1),path(R2)
output:
	tuple val(meta),path("OUT/*.gz"),emit:fastqs
	path("versions.yml"),emit:versions
script:
	def arsgs1=task.args.args1?:"--by-length 3G"
"""
set -o pipefail

mkdir -p TMP

seqkit split2 \\
	${args1} \\
	--force \\
	--extension ".gz" \\
	-O TMP --force \\
    --threads ${task.cpus} \\
    --read1 '${R1}' \\
    ${R2?"--read2 ${R2}":""}

mv -v TMP OUT

find OUT -type f 1>&2


touch versions.yml
"""
}
