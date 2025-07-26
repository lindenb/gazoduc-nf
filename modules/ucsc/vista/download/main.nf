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

process VISTA_DOWNLOAD {
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:bed
	tuple val(meta1),path("*.md"),emit:doc
	tuple path("versions.yml"),emit:versions
script:
   	def TAG = "VISTA"
	def whatis="VISTA enhancers"
	if(!meta1.ucsc_name || meta1.ucsc_name.isEmpty()) throw new IllegalArgumentException("${task.process} missing ucsc name");
	def url = "https://hgdownload.cse.ucsc.edu/gbdb/${meta1.ucsc_name}/vistaEnhancers/vistaEnhancers.bb"
"""
mkdir -p TMP/CACHE


curl -L -o TMP/jeter.bb "${url}"

bigBedToBed -udcDir=TMP/CACHE TMP/jeter.bb stdout |\\
        cut -f1,2,3,4 |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R ${fasta} --column 1 --convert SKIP  |\\
        LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
        bgzip > TMP/${TAG}.bed.gz

tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}">' > ${TAG}.header

cat << EOF > ${TAG}.md
Vista enhancer.
EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}
