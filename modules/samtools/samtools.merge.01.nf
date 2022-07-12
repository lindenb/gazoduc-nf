/*

Copyright (c) 2022 Pierre Lindenbaum

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

include {parseBoolean;assertNotEmpty;getKeyValue;moduleLoad} from '../../modules/utils/functions.nf'

process SAMTOOLS_MERGE_01 {
tag "${sample} N=${L.size()}"
afterScript "rm -rf TMP"
memory "5g"
cpus 5
input:
	val(meta)
	val(reference)
	tuple val(sample),val(L)
output:
	tuple val(sample),path("${sample}.merged.bam"),emit:bam
	path("version.xml"),emit:version
script:
	def with_index = parseBoolean(meta.with_index)
	def contig = (meta.contig?" -R \"${meta.contig}\" ":"")
	def interval = (meta.interval?" -R \"${meta.interval}\" ":"")
"""
hostname 2>&1
${moduleLoad("samtools")}
set -o pipefail
mkdir TMP

samtools merge --reference "${reference}" ${contig} ${interval} --threads ${task.cpus} -o TMP/jeter.bam ${L.join(" ")}

if [ ! -z "${with_index?"Y":""}"  ] ; then
	samtools index -@ ${task.cpus} TMP/jeter.bam
fi

mv TMP/jeter.bam "${sample}.merged.bam"

if [ ! -z "${with_index?"Y":""}"  ] ; then
	mv TMP/jeter.bam.bai "${sample}.merged.bam.bai" 
fi

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to fastq</entry>
	<entry key="sample">${sample}</entry>
	<entry key="count(bam)">${L.size()}</entry>
        <entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
</properties>
EOF
"""
stub:
"""
touch "${sample}.merged.bam" "${sample}.merged.bam.bai"
echo "<properties/>" > version.xml
"""
}

