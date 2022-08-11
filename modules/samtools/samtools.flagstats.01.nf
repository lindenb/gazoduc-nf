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

String getExtension(String fmt) {
	if(fmt==null) return ".txt";
	if(fmt.equals("tsv")) return ".tsv";
	if(fmt.equals("json")) return ".json";
	return ".txt";
	}

process SAMTOOLS_FLAGSTATS_01 {
tag "${row.bam}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${row.prefix?:""}flagstats${getExtension(meta.format)}"),emit:output
	path("version.xml"),emit:version
script:
	def prefix = row.prefix?:""
	def bed = row.bed==null || row.bed.name.equals("NO_FILE")?"":"-M -L \"${row.bed}\" "
	def interval = row.interval?:""
	def ref = row.reference?"--reference \"${row.reference}\" ":""
	def has_interval = !isBlank(bed) || !isBlank(interval)

"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

${has_interval?"samtools view --uncompressed -O BAM ${bed} ${ref} \"${row.bam}\" ${interval} \\":""}
samtools flagstats --threads ${task.cpus} -O "${row.format?:"default"}" ${has_interval?"-":"\"${row.bam}\""}  > "${prefix}flagstats${getExtension(meta.format)}"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">samtools flagtats</entry>
	<entry key="bam">${row.bam}</entry>
	<entry key="samtools.version">\$(samtools version | head -n1 )</entry>
</properties>
EOF
"""
}
