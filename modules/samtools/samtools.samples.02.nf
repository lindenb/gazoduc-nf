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

include { getKeyValue; moduleLoad; parseBoolean} from '../../modules/utils/functions.nf'

process SAMTOOLS_SAMPLES02 {
tag "${row.bams.name}"
afterScript "rm -f jeter.txt jeter.tsv"
input:
	val(meta)
	val(row)
output:
	path("sample2bam.tsv"),emit:output
	path("version.xml"),emit:version
script:
	def references = row.references==null || row.references.name.equals("NO_FILE")?"":" -F \"${row.references}\" "
	def with_header = parseBoolean(meta.with_header?:"true")
"""
hostname 2>&1
${moduleLoad("samtools")}
set -o pipefail

samtools samples \
	-f "${row.reference}" \
	${references} < "${row.bams}" |\
	sort -T . -t '\t' -k1,1 | uniq |\
	awk -F '\t' 'function uniq(S) {P=S;i=1; while(P in U) {i++;P=sprintf("%s.%d",S,i);} U[P]++; return P;} {printf("%s\t%s\t%s\t%s\\n",\$1,uniq(\$1),\$2,\$3);}' > jeter.tsv

# no empty samples
awk -F '\t' '(\$1==".")' jeter.tsv > jeter.txt 
test ! -s jeter.txt

# no empty ref
awk -F '\t' '(\$4==".")' jeter.tsv > jeter.txt 
test ! -s jeter.txt

# no dup samples
cut -f 2 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt

# no dup bam
cut -f 3 jeter.tsv | sort -T . | uniq -d > jeter.txt
test ! -s jeter.txt

awk -F '\t'  '(\$1!=\$2)' jeter.tsv > duplicate.names.txt


if [ ! -z "${with_header?"Y":""}" ] ; then
	echo "sample\tnew_sample\tbam\treference" > sample2bam.tsv
fi

cat jeter.tsv >>  sample2bam.tsv


##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract sample names from BAM metadata</entry>
	<entry key="input">${row.bams}</entry>
        <entry key="samtools.version">\$(samtools  --version | head -n 1| cut -d ' ' -f2)</entry>
        <entry key="n-samples">\$(wc -l < sample2bam.tsv )</entry>
        <entry key="samples">\$(cut -f 1 sample2bam.tsv |paste -s -d ' ')</entry>
        <entry key="duplicate names"><pre>\$(cut -f1,3 duplicate.names.txt)</pre></entry>
        <entry key="distinct references"><pre>\$(cut -f4 sample2bam.tsv | sort | uniq)</pre></entry>
</properties>
EOF
"""
stub:
"""
echo -e "sample\tnew_sample\tbam\treference" >  sample2bam.tsv
echo "<properties/>" > version.xml
"""
}
