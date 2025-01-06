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

include {getVersionCmd;isHg19;isHg38;moduleLoad} from '../utils/functions.nf'


process DOWNLOAD_REFGENE {
tag "${fasta}"
afterScript "rm -rf TMP"
input:
	path(fasta)
	path(fai)
	path(dict)
output:
	path("${fasta.baseName}.refGene.txt.*"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit htslib")}
mkdir -p TMP
set -o pipefail


cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:248956422	https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
1:249250621	https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
1:123556469	https://hgdownload.cse.ucsc.edu/goldenPath/canFam4/database/refGene.txt.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url


wget -O TMP/refGene.txt.gz `cat TMP/jeter.url`
gunzip -c TMP/refGene.txt.gz |\
		java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${fasta}" --column 3 --convert SKIP |\
		LC_ALL=C sort -T TMP -t '\t' -k3,3 -k5,5n |\
		bgzip > "${fasta.baseName}.refGene.txt.gz"

tabix -f -0 -b 5 -e 6 -s 3 "${fasta.baseName}.refGene.txt.gz"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download UCSC refgene</entry>
	<entry key="versions">${getVersionCmd("jvarkit/bedrenamechr wget")}</entry>
</properties>
EOF
"""
}

