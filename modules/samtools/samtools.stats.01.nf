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
if(!params.containsKey("gzip")) throw new IllegalArgumentException("params.gzip is missing");

include {parseBoolean;moduleLoad;isBlank} from '../../modules/utils/functions.nf'



process SAMTOOLS_STATS_01 {
tag "${row.sample}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(row)
output:
	tuple val(row),path("*.stats.txt*"),emit:output
	path("version.xml"),emit:version
script:
	def prefix = row.prefix?:(params.prefix?:"")
	def extra = row.extraSamtoolsStats?:"--remove-dups "
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

samtools stats ${extra}  --ref-seq "${row.reference}" --reference "${row.reference}" "${row.bam}" |\
	awk '/^# The command line/ {printf("# sample : ${row.sample}\\n");} {print}' > "${prefix}${row.sample}.stats.txt"

if ${parseBoolean(params.gzip)} ; then
gzip --best "${prefix}${row.sample}.stats.txt"
fi

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">samtools stats</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="bam">${row.bam}</entry>
	<entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}
