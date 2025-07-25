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
include {k1_signature} from '../../../modules/utils/k1.nf'

process PANMASK_DOWNLOAD{
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:output
	path("versions.yml"),emit:versions
	path("doc.md"),emit:doc
script:
	def k1 = k1_signature()
    def base = "https://zenodo.org/records/15683328/files/"
	def TAG = task.ext.tag?:"PANMASK"

"""
hostname 1>&2
mkdir -p TMP

export LC_ALL=C
set -x

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\t${base}/hg38.pm151a-v2.easy.bed.gz?download=1
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sed 's/^chr//' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
	sort | uniq > TMP/jeter.url



URL=`cat TMP/jeter.url`

if test -s TMP/jeter.url
then

	cut -f1,2 ${fai} | sort  -T TMP -t '\t' -k1,1 -k2,2n > TMP/genome.size

	wget --no-check-certificate -O - "\${URL}" |\\
		gunzip -c |\\
		cut -f1,2,3 |\\
		jvarkit -XX:-UsePerfData  -Djava.io.tmpdir=TMP bedrenamechr -f "${fasta}" --column 1 --convert SKIP  |\\
		sort  -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools complement -i - -g TMP/genome.size |\\
		awk -F '\t' '(\$1 ~ /^(chr)?[0-9XY]+\$/)' |\\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\\
		sed 's/\$/\t1/' |\\
		bgzip  > TMP/${TAG}.bed.gz

	echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="Panmask provides a list of hard regions for short-read variant calling against the human genome GRCh38 . See https://zenodo.org/records/15683328">' > ${TAG}.header


else
	touch TMP/${TAG}.bed
	bgzip TMP/${TAG}.bed

	echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="Panmask was NOT available for this build">' > ${TAG}.header

fi

tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./


cat << 'EOF' > doc.md
# annotations:panmask

> Panmask (https://zenodo.org/records/15683328) provides a list of easy/hard regions for short-read variant calling against the 
> human genome GRCh38. 

The flag ` TAG` means that the region is hard to call.

EOF


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "\${URL}"
END_VERSIONS
"""
}

