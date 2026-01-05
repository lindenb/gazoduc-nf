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

process TISSUES_DOWNLOAD {
tag "${meta1.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(gtf)
output:
    tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"), emit:bed
    path("versions.yml"),emit:versions
    path("doc.md"),emit:doc
script:
    def TAG = task.ext.tag?:"TISSUES"
    def URL = task.ext.url?:"https://download.jensenlab.org/human_tissue_knowledge_full.tsv"
    def WHATIZ = "Tissue from https://tissues.jensenlab.org/Search  (${URL})"
    def prefix = task.ext.prefix?:"${meta1.id}.${TAG}"
"""
hostname 1>&2
mkdir -p TMP
export LC_ALL=C

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name"   "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	cut -f1,2,3,5 |\\
        LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed

wget -O - "${URL}" |\\
    cut -f2,3,4 |\\
    awk -F '\t' '{OFS="\t";gsub(/[^A-Za-z0-9_:]+/,"_",\$3);print}' |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 |\\
    uniq > TMP/genes.txt

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2,2.3' TMP/genes.bed TMP/genes.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${prefix}.bed.gz

tabix --force -p bed TMP/${prefix}.bed.gz


mv TMP/*.bed.gz ./
mv TMP/*.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ}">' > ${prefix}.header
echo '##INFO=<ID=${TAG}_BRENDA,Number=.,Type=String,Description="${WHATIZ}">' >> ${prefix}.header


cat << 'EOF' > doc.md
# annotations:tissues

> TISSUES (https://tissues.jensenlab.org/Search) is a weekly updated web resource that integrates evidence on tissue expression 
> from manually curated literature, proteomics and transcriptomics screens, and automatic text mining. 

Tags:

 - `INFO/${TAG}` is the name of the tissue
 - `INFO/${TAG}_BRENDA` is the identifier in the "Brenda Ontology".

> The BRENDA tissue ontology (BTO) represents a comprehensive structured encyclopedia.
> It provides terms, classifications, and definitions of tissues, organs, anatomical structures, 
> plant parts, cell cultures, cell types, and cell lines of organisms from all taxonomic groups
> (animals, plants, fungi, protozoon) as enzyme sources. 

EOF


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${URL}"
END_VERSIONS
"""

stub:
    def prefix = task.ext.prefix?:"tissues"
"""
touch versions.yml ${prefix}.header ${prefix}.bed.gz ${prefix}.bed.gz.tbi doc.md
"""
}
