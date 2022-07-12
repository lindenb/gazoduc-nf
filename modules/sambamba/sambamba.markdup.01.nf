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

include {parseBoolean;moduleLoad} from '../../modules/utils/functions.nf'

process SAMBAMBA_MARKDUP_01 {
tag "${sample} ${file(bam).name}"
afterScript "rm -rf TMP"
memory "5g"
cpus 5
input:
	val(meta)
	tuple val(sample),val(bam)
output:
	tuple val(sample),path("${sample}.markdup.bam"),emit:bam
	path("version.xml"),emit:version
script:
	def with_index = parseBoolean(meta.with_index?:"true")
"""
hostname 2>&1
${moduleLoad("sambamba/0.0.0 samtools")}
set -o pipefail
mkdir TMP

# markdup
sambamba markdup --tmpdir=TMP -t ${task.cpus} "${bam}" TMP/jeter2.bam
mv TMP/jeter2.bam TMP/jeter.bam

if [ ! -z "${with_index?"Y":""}" ] ; then
	samtools index  -@ ${task.cpus} "TMP/jeter.bam"  "TMP/jeter.bam.bai"
fi

mv TMP/jeter.bam "${sample}.markdup.bam"

if [ ! -z "${with_index?"Y":""}" ] ; then
	mv TMP/jeter.bam.bai "${sample}.markdup.bam.bai"
fi

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">mardup with sambamba</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${bam}</entry>
        <entry key="sambamba.version">\$(sambamba --version 2>&1 | sort | uniq | paste -s -d ' ')</entry>
        <entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
</properties>
EOF
"""
stub:
"""
touch  "${sample}.markdup.bam"   "${sample}.markdup.bam.bai"
echo "<properties/>" > version.xml
"""
}

