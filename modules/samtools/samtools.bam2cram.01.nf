/*

Copyright (c) 2023 Pierre Lindenbaum

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

include {moduleLoad} from '../../modules/utils/functions.nf'

process SAMTOOLS_BAM_TO_CRAM_01 {
tag "${sample}"
afterScript "rm -rf TMP"
cpus 5
input:
	val(meta)
	val(reference)
	tuple val(sample),val(bam)
output:
	tuple val(sample),path("${meta.prefix?:""}${sample}.cram"),emit:cram
	path("${meta.prefix?:""}${sample}.cram.crai"),emit:index
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools")}

samtools view -@ ${task.cpus} -O CRAM -o "${meta.prefix?:""}${sample}.cram" -T "${reference}" "${bam}"
samtools index -@ ${task.cpus}  "${meta.prefix?:""}${sample}.cram"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to CRAM</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${bam}</entry>
	<entry key="size.in">\$(ls -la "${bam}")</entry>
	<entry key="size.out">\$(ls -la *.cram)</entry>
        <entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
</properties>
EOF
"""
stub:
"""
touch "${meta.prefix?:""}${sample}.cram"
touch "${meta.prefix?:""}${sample}.cram.crai"
echo "<properties/>" > version.xml
"""
}

