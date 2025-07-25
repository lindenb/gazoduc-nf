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
process REVEL_DOWNLOAD {
tag "${meta1.id?:fasta.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:output
	path("versions.yml"),emit:versions
	path("doc.md"),emit:doc
script:
	def TAG = "REVEL"
	def WHATIZ="REVEL is an ensemble method for predicting the pathogenicity of missense variants in the human genome. https://doi.org/10.1016/j.ajhg.2016.08.016 ."
	def url = task.ext.url?:"https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip?download=1"
	def colpos= -1;
	if(meta1.ucsc_name && meta1.ucsc_name.equals("hg19")) {
		colpos=2;
		}
	else  if(meta1.ucsc_name && meta1.ucsc_name.equals("hg38")) {
		colpos=3;
		}
	else {
		throw new IllegalArgumentException("Bad fai version ${task.process}");
		}
"""
hostname 1>&2
mkdir -p TMP
wget -O TMP/jeter.zip "${url}"


unzip -p TMP/jeter.zip  revel_with_transcript_ids |\\
	tail -n +2 |\\
	tr "," "\t" |\\
	cut -f 1,${colpos},4,5,8 |\\
	awk -F '\t' '\$2!="."' |\\
	uniq |\\
	jvarkit  -Djava.io.tmpdir=TMP  bedrenamechr -R "${fasta}" --column 1 --convert SKIP  |\\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n |\\
		uniq |\\
		bgzip > TMP/${TAG}.bed.gz

tabix -s 1 -b 2 -e 2 -f TMP/${TAG}.bed.gz


mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=1,Type=Float,Description="${WHATIZ} ${url}">' > ${TAG}.header


cat << 'EOF' > doc.md
# annotations:revel

`INFO/${TAG}` : REVEL

> REVEL is an ensemble method for predicting the pathogenicity of missense variants in the 
> human genome. https://doi.org/10.1016/j.ajhg.2016.08.016

EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}
