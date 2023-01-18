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
include {parseBoolean;moduleLoad;getVersionCmd} from '../utils/functions.nf'

/**
 split one bed file into one bed per contig
*/
process BED_SPLIT_CONTIG_01 {
tag "${bed.name}"
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("contigs.bed"),emit:bed /* chrom / start / end / path to bed */
	path("contigs.txt"),emit:contigs /* unique contigs */
	path("version.xml"),emit:version
script:
	def with_header= parseBoolean(meta.with_header?:false)
	def with_tabix= parseBoolean(meta.with_tabix?:false)
"""
hostname 1>&2
${moduleLoad("htslib")}
set -o pipefail

mkdir -p BEDS

${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\
	grep -vE "^(#|browser|track)" |\
	sort -t '\t' -k1,1 -k2,2n -T . |\
	awk -F '\t' 'BEGIN{C="";N=0;} {OFS="\t";if(C!=\$1) {if(N>0) {close(f);} N=0;C=\$1;f=sprintf("BEDS/%s.bed",C);} print >> f; N++;} END {if(N>0) close(f);}'

touch jeter.bed



find \${PWD}/BEDS -type f -name "*.bed" | while read F
do
	awk -F '\t' -vF=\$F -vX="${with_tabix?".gz":""}" '{C=\$1;b=int(\$2);e=int(\$3);if(NR==1 || b<B) B=b; if(NR==1 || e>E) E=e;} END{printf("%s\t%d\t%d\t%s%s\\n",C,B,E,F,X);}' "\${F}" >> jeter.bed

	if ${with_tabix} ; then
		bzgip "\${F}"
		tabix -p bed "\${F}"
	fi
done


if ${with_header} ; then
	echo -e "contig\tstart\tend\tvcf" > contigs.bed
fi

sort -T . -t '\t' -k1,1 -k2,2n jeter.bed >> contigs.bed
rm jeter.bed

cut -f1 contigs.bed | sort | uniq > contigs.txt

#######################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">split bed file to one bed per contig</entry>
	<entry key="bed">${bed.toRealPath()}</entry>
	<entry key="version">${getVersionCmd("awk tabix")}</entry>
</properties>
EOF
"""
}
