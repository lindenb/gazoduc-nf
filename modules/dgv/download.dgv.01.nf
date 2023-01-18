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

include {isBlank;moduleLoad;isHg38;isHg19} from '../utils/functions.nf'

process DOWNLOAD_DGV_01 {
input:
	val(meta)
	val(reference)
output:
       	path("dgv.bed.gz"),emit:bed
       	path("dgv.bed.gz.tbi"),emit:index
        path("version.xml"),emit:version
script:
	def url=isHg19(reference)?"http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt":(isHg38(reference)?"http://dgv.tcag.ca/dgv/docs/GRCh38_hg38_variants_2020-02-25.txt":"")
	def whatis ="DGV (Database of Genomic Variants) as BED file from ${url}"
if(!isBlank(url))
"""
hostname 1>&2
${moduleLoad("htslib jvarkit bedtools")}
set -o pipefail
mkdir TMP
wget  -O - "${url}" |\
	awk -F '\t' '{printf("%s\t%s\t%s\t%s",\$2,\$3,\$4,\$1);for(i=5;i<=NF;i++){printf("\t%s",\$i);}printf("\\n");}' |\
	sed  's/^chr/#chr/' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/dgv.bed


bgzip TMP/dgv.bed
tabix --comment '#'  -p bed -f TMP/dgv.bed.gz

mv TMP/dgv.bed.gz ./
mv TMP/dgv.bed.gz.tbi ./

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
else
"""
${moduleLoad("htslib")}
touch dgv.bed
bgzip dgv.bed
tabix --comment '#' -p bed -f dgv.bed.gz
echo '<properties/>' > version.xml

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">${whatis}</entry>
        <entry key="url">no DGV is available for ${reference}</entry>
</properties>
EOF
"""
}
