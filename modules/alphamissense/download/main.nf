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

process ALPHAMISSENSE_DOWNLOAD {
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(dict)
output:
	tuple val(meta1),path("*.tsv.gz"),path("*.tsv.gz.tbi"),path("*.header"),emit:output
	tuple val(meta1),path("*.md"),emit:doc
	path("versions.yml"),emit:versions
script:
    def TAG = task.ext.tag?:"ALPHAMISSENSE"
    def whatis="Data from AlphaMissense https://www.science.org/doi/10.1126/science.adg7492 Accurate proteome-wide missense variant effect prediction with AlphaMissense"
    def base="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_"
    if(!meta1.ucsc_name) throw new IllegalArgumentException("${task.process} undefined ucsc_name");
    def url="";
    if(meta1.ucsc_name.equals("hg38")) {
        url = base+"hg38.tsv.gz";
    } else if(meta1.ucsc_name.equals("hg19")) {
        url = base+"hg19.tsv.gz";
    } else {
        url=""
        }
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

if ${url.isEmpty()}
then
	touch TMP/${TAG}.tsv
	bgzip TMP/${TAG}.tsv
else

	curl -L -o TMP/jeter.tsv.gz "${url}"

	gunzip -c TMP/jeter.tsv.gz  |\\
		grep -v '^#' |\\
		cut -f 1 | uniq | LC_ALL=C sort -T TMP | uniq |\\
		awk '{printf("%s\t%s\\n",\$1,\$1);}' |\\
		jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -f "${dict}" --column 2 --convert SKIP |\\
		awk -F '\t' '{printf("s|^%s\t|%s\t|\\n",\$1,\$2);}' > TMP/jeter.sed


	gunzip -c TMP/jeter.tsv.gz |\\
		grep -v '^#'  |\\
		cut -f 1-4,9,10 |\\
		sed -f TMP/jeter.sed |\\
		LC_ALL=C sort --buffer-size=${task.memory.mega}M -T TMP -t '\t' -k1,1 -k2,2n |\\
		uniq |\\
		bgzip >  TMP/${TAG}.tsv.gz 

fi

tabix -s 1 -b 2 -e 2  TMP/${TAG}.tsv.gz
mv TMP/${TAG}.tsv.gz ./
mv TMP/${TAG}.tsv.gz.tbi ./

echo '##INFO=<ID=${TAG}_PATHOGENOCITY,Number=1,Type=Float,Description="${whatis}.">' >  ${TAG}.header
echo '##INFO=<ID=${TAG}_CLASS,Number=.,Type=String,Description="${whatis}.">' >> ${TAG}.header


cat << EOF > ${TAG}.md
${whatis}
EOF

cat << EOF > versions.yml
${task.process}:
	jvarkit: TODO
	url : ${url}
EOF
"""

stub:
	def TAG = task.ext.tag?:"ALPHAMISSENSE" 
"""
touch versions.yml ${TAG}.header ${TAG}.tsv.gz ${TAG}.tsv.gz.tbi ${TAG}.md
"""
}
