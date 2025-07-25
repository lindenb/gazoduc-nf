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

process RMSK_DOWNLOAD{
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:bed
	path("versions.yml"),emit:versions
	path("doc.md"),emit:doc
script:
    if(!meta1.ucsc_name || meta1.ucsc_name.isEmpty()) throw new IllegalArgumentException("${task.process} no meta1.ucsc_name");
    def url = "https://hgdownload.cse.ucsc.edu/goldenPath/${meta1.ucsc_name}/database/rmsk.txt.gz"
    def TAG = task.ext.tag?:"RMSK"

"""
hostname 1>&2
mkdir -p TMP

set -o pipefail
curl -L  "${url}" |\\
	gunzip -c |\\
	cut -f6-8 |\\
	jvarkit bedrenamechr -XX:-UsePerfData  -Djava.io.tmpdir=TMP -f "${fasta}" --column 1 --convert SKIP  |\\
		sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge |\\
		sed 's/\$/\t1/' |\\
		bgzip > TMP/${TAG}.bed.gz
	

tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="Repeat Masker from UCSC">' > ${TAG}.header

cat << 'EOF' > doc.md
# annotations:repeatmasker

`INFO/${TAG}` : Repeat Masker regions from UCSC

EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}

