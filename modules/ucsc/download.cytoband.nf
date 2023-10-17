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

include {getVersionCmd;isHg19;isHg38;moduleLoad} from '../utils/functions.nf'


process DOWNLOAD_CYTOBAND {
tag "${reference}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
output:
	path("${genomeId}.cytoBandIdeo.txt"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url="https://hgdownload.cse.ucsc.edu/goldenPath/${genome.ucsc_name}/database/cytoBandIdeo.txt.gz"
"""
hostname 1>&2
${moduleLoad("jvarkit")}
mkdir -p TMP
set -o pipefail


wget -O TMP/cytoBandIdeo.txt.gz "${url}"
gunzip -c TMP/cytoBandIdeo.txt.gz |\
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP > "${genomeId}.cytoBandIdeo.txt" 

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download UCSC cytobands</entry>
	<entry key="reference">${reference}</entry>
	<entry key="url"><url>${url}</url></entry>
	<entry key="versions">${getVersionCmd("jvarkit/bedrenamechr wget")}</entry>
</properties>
EOF
"""
}

