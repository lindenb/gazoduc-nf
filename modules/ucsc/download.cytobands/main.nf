/*
Copyright (c) 2025 Pierre Lindenbaum

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


process DOWNLOAD_CYTOBAND {
tag "${meta1.id?:fasta.name}"
label "process_quick"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.cytoBandIdeo.txt"),emit:output
	path("versions.yml"),emit:versions
script:
	def k1 = k1_signature()
	def base = "https://hgdownload.cse.ucsc.edu/goldenPath"
	def prefix = task.ext.prefix?:fasta.baseName
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg19}\t${base}/hg19/database/cytoBandIdeo.txt.gz
1:${k1.hg38}\t${base}/hg38/database/cytoBandIdeo.txt.gz
1:${k1.canFam3}\t${base}/canFam3/database/cytoBandIdeo.txt.gz
1:${k1.canFam4}\t${base}/canFam4/database/cytoBandIdeo.txt.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sed 's/^chr//' |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
	sort | uniq > TMP/jeter.url

test -s TMP/jeter.url


wget -O TMP/cytoBandIdeo.txt.gz `cat TMP/jeter.url`

gunzip -c TMP/cytoBandIdeo.txt.gz |\\
		jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP > "${prefix}.cytoBandIdeo.txt" 

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "\${`cat TMP/jeter.url`}"
END_VERSIONS
"""
}

