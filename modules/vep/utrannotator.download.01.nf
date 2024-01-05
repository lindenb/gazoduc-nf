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
include {getVersionCmd;moduleLoad} from '../utils/functions.nf'
/**

A VEP Plugin to annotate high-impact five prime UTR variants either creating new upstream ORFs or disrupting existing upstream ORFs

*/
process UTR_ANNOTATOR_DOWNLOAD_01 {
tag "${genomeId}"
input:
	val(meta)
	val(genomeId)
output:
	path("UTRannotator/uORF_5UTR.txt"),emit:output
	path("UTRannotator/uORF_5UTR.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome  = params.genomes[genomeId]
	def reference = genome.fasta
	def slop = 20
	def url = "https://github.com/ImperialCardioGenetics/UTRannotator"
	def build = genome.ensembl_name
"""
hostname 1>&2
set -o pipefail
${moduleLoad("bedtools  jvarkit")}

test ! -z "${build}"

git clone "${url}.git"

ln -s "\${PWD}/UTRannotator/uORF_5UTR_${build}_PUBLIC.txt" "\${PWD}/UTRannotator/uORF_5UTR.txt"

awk -F '\t' '(\$1!="chr") {B0=int(\$2);E0=int(\$6);if(B0<E0) {B=B0;E=E0;} else {B=E0;E=B0;} printf("%s\t%d\t%d\\n",\$1,B,E);}' "\${PWD}/UTRannotator/uORF_5UTR.txt" |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP |\
	bedtools slop -b ${slop} -g "${reference}.fai" |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > UTRannotator/uORF_5UTR.bed

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download VEP plugin for UTRs</entry>
	<entry key="url"><url>${url}</url></entry>
	<entry key="reference">${reference}</entry>
	<entry key="slop">${slop}</entry>
	<entry key="version">${getVersionCmd("git bedtools")}</entry>
</properties>
EOF
"""
}
