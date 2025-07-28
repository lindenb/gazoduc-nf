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

process ENSEMBL_REG_DOWNLOAD{
tag "${fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:bed
	path("versions.yml"),emit:versions
script:
   	def TAG = "ENSEMBL_REG"
	def whatis="Ensembl regulatory Features"
    f(!meta1.ucsc_name) throw new IllegalArgumentException("${task.process} undefined ucsc_name");
    def url="";
    if(meta1.ucsc_name.equals("hg38")) {
        url = "https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v114.gff3.gz"
    } else if(meta1.ucsc_name.equals("hg19")) {
        url = "https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh37/annotation/Homo_sapiens.GRCh37.regulatory_features.v114.gff3.gz"
    } else {
        throw new IllegalArgumentException("${task.process} unknown ucsc_name");
        }
"""
hostname 1>&2
mkdir -p TMP

curl -L "${url}" |\\
        gunzip -c |\\
	    awk -F '\t' '/^[^#]/{printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\\
        jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R ${fasta} --column 1 --convert SKIP  |\\
        LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	    uniq | bgzip > TMP/${TAG}.bed.gz
	
tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis} ${url}">' > ${TAG}.header

cat << EOF > versions.yml
${task.process}:
	url: "${url}"
EOF
"""
}