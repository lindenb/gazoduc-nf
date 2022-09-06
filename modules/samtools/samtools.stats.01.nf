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

include {assertNotEmpty;getKeyValue;moduleLoad;isBlank} from '../../modules/utils/functions.nf'

process SAMTOOLS_STATS_01 {
tag "${row.sample}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${row.sample}.stats.tsv.gz"),emit:output
	path("version.xml"),emit:version
script:
	def  extra = ${meta.extraSamtoolsStats?:""}
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

samtools stats ${extra}  --ref-seq "${row.reference}" --reference "${row.reference}" "${row.bam}" |\
	awk '/^# The command line/ {printf("# sample : ${row.new_sample}\\n");} {print}' > "${row.sample}.stats"

gzip --best "${row.sample}.stats"

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
