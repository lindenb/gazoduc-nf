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


/*
Example output:

S1.R1.part_001.fq.gz    S1.R1.part_002.fq.gz
S1.R1.part_003.fq.gz    S1.R1.part_004.fq.gz
S1.R1.part_005.fq.gz    S1.R2.part_001.fq.gz
S1.R2.part_002.fq.gz    S1.R2.part_003.fq.gz
S1.R2.part_004.fq.gz    S1.R2.part_005.fq.gz


*/
process SEQKIT_SPLIT2 {
label "process_single"
conda "${moduleDir}/../../../conda/seqkit.yml"
tag "${meta.id} ${R1.name} ${R2?R2.name:""}"
afterScript "rm -rf TMP"
// cpus 4 set in config file
input:
	tuple val(meta),path(R1),path(R2)
output:
	tuple val(meta),path("OUT/*.gz",arity: '1..*'),emit:fastqs
	path("versions.yml"),emit:versions
script:
	def args1 = (task.ext.args1?:"")
	if(args1.toString().isEmpty()) throw new IllegalArgumentException("${task.process} undefined task.args1 e.g. --by-length 3G");
"""
mkdir -p TMP

seqkit split2 \\
	${args1} \\
	--force \\
	--extension ".gz" \\
	-O TMP \\
	--force \\
    --threads ${task.cpus} \\
    --read1 "${R1}" \\
    ${R2?"--read2 \"${R2}\"":""}

mv -v TMP OUT

find OUT -type f 1>&2


cat <<-END_VERSIONS > versions.yml
"${task.process}":
        seqkit: todo
END_VERSIONS
"""
stub:
	def args1 = (task.ext.args1?:"")
	if(args1.toString().isEmpty()) throw new IllegalArgumentException("${task.process} undefined task.args1 e.g. --by-length 3G");
"""
touch ${R1.name}.part_001.fq ${R1.name}.part_002.fq ${R2?"{R2.name}.part_001.fq ${R2.name}.part_002.fq":""}
gzip -f *.fq
touch versions.yml 
"""
}
