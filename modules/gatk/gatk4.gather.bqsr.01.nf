/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {moduleLoad} from './../utils/functions.nf'

process GATK4_GATHER_BQSR_01 {
tag "${key.sample} ${key.bam} N=${L.size()}"
cache 'lenient'
memory '5g'
cpus 1
afterScript 'rm -rf TMP'
input:
	val(meta)
	tuple val(key),val(L)
output:
        tuple val(key),path("${key.sample}.recal.table"),emit:output
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("gatk4")}

mkdir -p TMP

# allele specific annotation are not supported in non-gvcf mode
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GatherBQSRReports \
	${L.collect(V->"-I "+V).join(" ")} \
	-O "${key.sample}.recal.table" 

###########################################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">gatk GatherBQSRReports</entry>
	<entry key="sample">${key.sample}</entry>
	<entry key="bam">${key.bam}</entry>
	<entry key="count">${L.size()}</entry>
	<entry key="gatk.version">\$( gatk --version 2> /dev/null  | paste -s -d ' ' )</entry>
</properties>
EOF
"""

stub:
"""
touch "${key.sample}.recal.table"
echo "<properties/>" > version.xml
"""
}
