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

process RMSK_DOWNLOAD{
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:bed
	path("versions.yml"),emit:versions
	path("doc.md"),emit:doc
script:
    def url =task.ext.url?:""
    if(url.isEmpty()) {
		if(meta.ucsc_name==null || meta.ucsc_name.isEmpty()) throw new IllegalArgumentException("${task.process} no meta1.ucsc_name");
		url="https://hgdownload.cse.ucsc.edu/goldenPath/${meta.ucsc_name}/database/rmsk.txt.gz"
		}
    def TAG = task.ext.tag?:"RMSK"
    def prefix = task.ext.prefix?:"${meta.id}.${TAG}"
"""
hostname 1>&2
mkdir -p TMP

set -o pipefail
curl -L  "${url}" |\\
	gunzip -c |\\
	cut -f6-8 |\\
	jvarkit bedrenamechr -XX:-UsePerfData  -Djava.io.tmpdir=TMP -f "${dict}" --column 1 --convert SKIP  |\\
		sort  -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		bedtools merge |\\
		sed 's/\$/\t1/' |\\
		bgzip > TMP/${prefix}.bed.gz
	

tabix -p bed -f TMP/${prefix}.bed.gz

mv TMP/${prefix}.bed.gz ./
mv TMP/${prefix}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="Repeat Masker from UCSC">' > ${prefix}.header

cat << 'EOF' > doc.md
# annotations:repeatmasker

`INFO/${TAG}` : Repeat Masker regions from UCSC

EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
 def prefix="rmsk"
"""
touch versions.yml doc.md "${prefix}.bed.gz" "${prefix}.bed.gz.tbi" ${prefix}.header
"""
}

