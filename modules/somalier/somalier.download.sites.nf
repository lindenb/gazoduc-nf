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
include {moduleLoad} from '../utils/functions.nf'
include {k1_signature} from '../utils/k1.nf'

def k1 = k1_signature()

process SOMALIER_DOWNLOAD_SITES {
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	path(genome)
output:
	path("sites.vcf.*"),emit:output
script:
	def fasta = genome.find{it.name.endsWith("a")}
	def fai = genome.find{it.name.endsWith(".fai")}
"""
hostname 1>&2
set -o pipefail
set -x
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg19}\t
EOF

cut -f1,2 "${fai}" |tr "\t" ":" | sort -T TMP  > TMP/jeter2.tsv

URL=`join -t '\t' -1 1 -2 1 -o "1.2" TMP/jeter1.tsv TMP/jeter2.tsv`

test ! -z "\${URL}"

wget -O - "\${URL}" |\\
	gunzip -c |\\
	jvarkit  vcfsetdict -R "${fasta}"  --onNotFound SKIP |\\
	bcftools sort -T TMP/sort -o TMP/sites.vcf.gz -O z

bcftools index -f -t TMP/sites.vcf.gz

mv -v TMP/sites* ./
"""
}
