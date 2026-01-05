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

include {k1_signature} from '../../utils/k1.nf'

process DOWNLOAD_GENCODE {
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	path(genome)
output:
	path("*.gencode.*"),emit:output
script:
	def k1 = k1_signature()
	def fai= genome.find{it.name.endsWith(".fai")}
	def dict= genome.find{it.name.endsWith(".dict")}
	def fasta= genome.find{it.name.endsWith("a")}
	def version="V47"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail


cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}	https://hgdownload.soe.ucsc.edu/goldenpath/hg38/database/wgEncodeGencodeComp${version}.txt.gz
1:${k1.hg19}	https://hgdownload.soe.ucsc.edu/goldenpath/hg19/database/wgEncodeGencodeComp${version}lift37.txt.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' | sed 's/^chr//' | sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv | sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

wget -O TMP/gencode.sql `cat TMP/jeter.url| sed 's/\\.txt\\.gz\$/.sql/'`

wget -O TMP/gencode.txt.gz `cat TMP/jeter.url`

gunzip -c TMP/gencode.txt.gz |\\
		jvarkit bedrenamechr -f "${fasta}" --column 3 --convert SKIP |\\
		LC_ALL=C sort -T TMP -t '\t' -k3,3 -k5,5n |\\
		bgzip > "${fasta.baseName}.gencode.txt.gz"

tabix -f -0 -b 5 -e 6 -s 3 "${fasta.baseName}.gencode.txt.gz"
mv TMP/gencode.sql  "${fasta.baseName}.gencode.sql"
"""
}

