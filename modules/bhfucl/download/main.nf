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

process BHFUCL_DOWNLOAD {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(gtf)
output:
	tuple val(meta),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"), emit:bed
	path("versions.yml"),emit:versions
script:
    def TAG = task.ext.tag?:"BHFUCL"
    def url = task.ext.url?:"http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
    def WHATIZ = "Cardiovascular Gene Ontology Annotation Initiative ${url}"
	def extend = task.ext.extend?:1000
    def prefix = task.ext.prefix?:"${meta.id}.${TAG}"
"""
hostname 1>&2
mkdir -p TMP
set -o pipefail
mkdir -p TMP

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name"  "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	cut -f1,2,3,5 |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed

curl -L "${url}" |\\
	gunzip -c |\\
	awk -F '\t' '/#/{next} (\$4 ~ /^NOT/) {next;} {print \$3;}' |\\
	uniq | LC_ALL=C sort -T TMP | uniq > TMP/genes.txt

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.4' TMP/genes.bed TMP/genes.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}.bed.gz

tabix --force -p bed TMP/${TAG}.bed.gz


mv TMP/*.bed.gz ./
mv TMP/*.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ} ${url}">' > ${prefix}.header


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml
touch ${prefix}.header ${prefix}.bed.gz ${prefix}.bed.gz.tbi
"""
}

