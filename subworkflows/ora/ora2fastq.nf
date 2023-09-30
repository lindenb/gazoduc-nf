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

/** convert ORA to FASTQ using illumina ORAD */
workflow ORA_TO_FASTQ {
	take:
		meta
		rows /** contains sample,R1 */
	main:
		
		version_ch = Channel.empty()
	
		rows.branch {
			fastq : !it.R1.endsWith(".ora")
			ora : it.R1.endsWith(".ora")
			}.set{branch_ch}
		

		ch1 = APPLY_ORAD([:] , branch_ch.ora)
		version_ch = version_ch.mix(ch1.version)

	
		ora_ch = branch_ch.ora.combine( ch1.output.splitCsv(sep:'\t',header:true) ).
			filter{T->T[0].R1.equals(T[1].ora)}.
			map{T->T[0].plus(T[1])}
	

		rows_out = branch_ch.fastq.mix(ora_ch)
	emit:
		version = version_ch
		output = rows_out
	}



process APPLY_ORAD {
tag "${row.sample} ${row.R1}"
afterScript "rm -rf TMP OUT/fastq.tsv"
cpus 4
input:
	val(meta)
	val(row)
output:
	path("fastq.tsv"),emit:output
	path("version.xml"),emit:version
script:
	if(row.containsKey("R2") && !row.R2.isEmpty() && !row.R2.equals(".")) throw new IllegalArgumentException("ORA specified but got R2 : ${row}");
"""
hostname 1>&2
${moduleLoad("orad")}

mkdir -p TMP
orad '${row.R1}' -V -N -t ${task.cpus} -P TMP
mv -v TMP OUT


find \${PWD}/OUT/ -type f -name "*q.gz" |\
	sort -T . -k1,1V |\
	paste - - |\
	awk -F '\t' '{printf("${row.R1}\t%s\\n",\$0);}' > OUT/fastq.tsv

test -s OUT/fastq.tsv

echo "ora\tR1\tR2" > fastq.tsv
cat OUT/fastq.tsv >> fastq.tsv



##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">orad2fastq</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="R1">${row.R1}</entry>
        <entry key="versions">\$(orad  --version)</entry>
</properties>
EOF
"""
}

