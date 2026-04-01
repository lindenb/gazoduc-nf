/*

Copyright (c) 2026 Pierre Lindenbaum

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

include {assertNotEmpty;getKeyValue;moduleLoad} from '../../modules/utils/functions.nf'

process SAMTOOLS_IDXSTATS_01 {
tag "{row.bam}"
afterScript "rm -rf TMP"
cpus 1
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("*.idxstats.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def prefix = row.prefix?:""
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

SN=`samtools samples "${row.bam}" | cut -f 1 | head -n 1`

test ! -z "${SN}"

samtools idxstats --threads ${task.cpus} "${row.bam}"  |\
	awk -F '\t' -vS=\${SN} 'BEGIN{OFS="\t";printf("sample\tbam\tcontig\tlength\tmapped\tunmapped\\n");}{printf("%s\t${bam}\t%s\\n",S,\$0);}' > "${prefix}\${SN}.idxstats.tsv"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">samtools idxstats</entry>
	<entry key="bam">${bam}</entry>
	<entry key="sample">\${SN}</entry>
	<entry key="mapped">\$(tail -n+2 "${prefix}\${SN}.idxstats.tsv" | cut -f 5  | paste -s -d '+' |bc)</entry>
	<entry key="unmapped">\$(tail -n+2 "${prefix}\${SN}.idxstats.tsv" | cut -f 6  | paste -s -d '+' |bc)</entry>
</properties>
EOF
"""
stub:
"""
touch  idxstats.tsv
echo "<properties/>" > version.xml
"""
}

