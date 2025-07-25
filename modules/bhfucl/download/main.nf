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

process BHFUCL_DOWNLOAD {
tag "${fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(gtf),path(gtf_tbi)
output:
	tuple val(meta1),path("*IN.bed.gz"), path("*IN.bed.gz.tbi"), path("*IN.header"), emit:bed
	tuple val(meta1),path("*OUT.bed.gz"),path("*OUT.bed.gz.tbi"),path("*OUT.header"),emit:bed_extended
	path("versions.yml"),emit:versions
script:
    def TAG = task.ext.tag?:"BHFUCL"
    def URL = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
    def WHATIZ = "Cardiovascular Gene Ontology Annotation Initiative ${URL}"
	def extend = task.ext.extend?:1000
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	cut -f1,2,3,5 |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed

wget -O - "${URL}" |\\
	gunzip -c |\\
	awk -F '\t' '/#/{next} (\$4 ~ /^NOT/) {next;} {print \$3;}' |\\
	uniq | LC_ALL=C sort -T TMP | uniq > TMP/genes.txt

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4' TMP/genes.bed TMP/genes.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}_IN.bed.gz

tabix --force -p bed TMP/${TAG}_IN.bed.gz

gunzip -c TMP/${TAG}_IN.bed.gz |\\
	bedtools slop -b ${extend} -g ${fai} |\\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP > TMP/extended.bed


bedtools subtract -a  TMP/extended.bed -b TMP/${TAG}_IN.bed.gz |\\
	awk -F '\t' 'int(\$2) < int(\$3)' |\\
	LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
	uniq |\\
	bgzip  > TMP/${TAG}_OUT.bed.gz

tabix --force -p bed TMP/${TAG}_OUT.bed.gz


mv TMP/*.bed.gz ./
mv TMP/*.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ} ${URL}">' > ${TAG}_IN.header
echo '##INFO=<ID=${TAG}_NEAR,Number=.,Type=String,Description="Near gene distance=${extend}. ${WHATIZ} ${URL}">' > ${TAG}_OUT.header


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${URL}"
END_VERSIONS
"""
}


