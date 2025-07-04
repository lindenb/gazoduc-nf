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

process DOWNLOAD_GTF_OR_GFF3 {
tag "${meta1.id?:fasta.name}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.gz"), path("*.gz.tbi"),emit:output
	path("versions.yml"),emit:versions
script:
	def k1 = k1_signature()
	def extension = (task.ext==null || task.ext.suffix==null?
		(task.process.toString().toLowerCase().endsWith("gtf")?"gtf":
			(task.process.toString().toLowerCase().endsWith("gff3")?"gff3":""))
		:task.ext.suffix
		)
	
	if(extension.isEmpty()) throw new IllegalArgumentException("suffix missing for ${task.process}");
	def release = task.ext.release?:"114"
	def base = task.ext.base?:"https://ftp.ensembl.org/pub"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail

##  GTF/GFF3 from gencode are not compatible with bcftools csq
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/release-${release}/${extension}/homo_sapiens/Homo_sapiens.GRCh38.${release}.${extension}.gz
1:${k1.hg19}\t${base}/grch37/current/${extension}/homo_sapiens/Homo_sapiens.GRCh37.87.chr.${extension}.gz
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sed 's/^chr//' |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
	sort |\\
	uniq > TMP/jeter.url

test -s TMP/jeter.url


wget -O TMP/gencode.txt.gz `cat TMP/jeter.url`

gunzip -c TMP/gencode.txt.gz |\\
		jvarkit bedrenamechr -f "${fasta}" --column 1 --convert SKIP |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k4,4n |\\
		bgzip > "${fasta.baseName}.${extension}.gz"

tabix -f -p gff  "${fasta.baseName}.${extension}.gz"

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(tabix version | awk '(NR==1) {print \$NF;}')"
	URL: "\$(cat TMP/jeter.url)"
END_VERSIONS
"""
}

