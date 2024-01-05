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
include {moduleLoad} from '../../modules/utils/functions.nf'

/** split FASTQ using seqtk split2 */
workflow SEQTK_SPLITFASTQ {
	take:
		meta
		rows /** contains sample,R1,(R2)? */
	main:
		if(!params.seqtk.split2.containsKey("args")) throw new IllegalArgumentException("params.seqtk.split2.args missing")
		
		version_ch = Channel.empty()
		
		

		if(params.seqtk.split2.args.trim().isEmpty()) {
			rows_out = rows
			}
		else
			{
			ch1 = SEQTK_SPLIT([:] ,rows)
			version_ch = version_ch.mix(ch1.version)

			rows_out = ch1.output.splitCsv(header:true, sep: '\t')
			}
		
		
	emit:
		version = version_ch
		output = rows_out
	}



process SEQTK_SPLIT {
tag "${row.sample} ${row.R1} ${row.R2?:""}"
afterScript "rm -rf TMP OUT/fastq.tsv"
// cpus 4 set in config file
input:
	val(meta)
	val(row)
output:
	path("fastq.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def COL = row.grep({!(it.key.equals("R1") || it.key.equals("R2"))}).collect{it.key}.join("\t")
	def ROW = row.grep({!(it.key.equals("R1") || it.key.equals("R2"))}).collect{it.value}.join("\t")

	seqtk = params.seqtk.split2.args

if(row.containsKey("R2") && !row.R2.isEmpty() && !row.R2.equals("."))
"""
hostname 1>&2
${moduleLoad("seqkit")}
set -o pipefail

mkdir -p TMP

seqkit split2 -O TMP --force -j ${task.cpus} ${seqtk} -1 '${row.R1}' -2 '${row.R2}'
mv -v TMP OUT


find \${PWD}/OUT/ -type f -name "*q.gz" |\
	awk -F '.' '{printf("%s\t%s\\n",\$(NF-2),\$0);}' |\
	sort -T . -t '\t' -k1,1V |\
	paste - - |\
	cut -f 2,4 |\
	awk -F '\t' '{printf("${ROW}\t%s\\n",\$0);}' > OUT/fastq.tsv

test -s OUT/fastq.tsv

echo "${COL}\tR1\tR2" > fastq.tsv
cat OUT/fastq.tsv >> fastq.tsv



##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">split fastq</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="R1">${row.R1}</entry>
	<entry key="R2">${row.R2}</entry>
	<entry key="seqtk">${seqtk}</entry>
        <entry key="versions">\$(seqkit  version)</entry>
</properties>
EOF
"""

else

"""
hostname 1>&2
${moduleLoad("seqkit")}
set -o pipefail

mkdir -p TMP

seqkit split2 -O TMP --force -j ${task.cpus} ${seqtk} ${row.R1}
mv -v TMP OUT


find \${PWD}/OUT/ -type f -name "*q.gz" |\
	sort -T . |\
	awk  '{printf("${ROW}\t%s\\n",\$0);}' > OUT/fastq.tsv

test -s OUT/fastq.tsv

echo "${COL}\tR1" > fastq.tsv
cat OUT/fastq.tsv >> fastq.tsv

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">split fastq</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="R1">${row.R1}</entry>
	<entry key="seqtk">${seqtk}</entry>
        <entry key="versions">\$(seqkit  version)</entry>
</properties>
EOF
"""
}

