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



workflow TISSUES {
	take:
		meta
		fasta
		fai
		dict
		gtf
		vcfs /* meta, vcf,vcf_index */
	main:
		versions = Channel.empty()

		DOWNLOAD(fasta,fai,dict,gtf)
		versions = versions.mix(DOWNLOAD.out.versions)

		ANNOTATE(
			DOWNLOAD.out.bed,
			vcfs
			)
		versions = versions.mix(ANNOTATE.out.versions)
	emit:
		
        vcf = ANNOTATE.out.vcf	
        versions
        doc = DOWNLOAD.out.doc

	}
		
process DOWNLOAD {
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
    tuple val(meta1),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"), emit:bed
    path("versions.yml"),emit:versions
    path("doc.md"),emit:doc
script:
    def TAG = task.ext.tag?:"TISSUES"
    def URL = task.ext.url?:"https://download.jensenlab.org/human_tissue_knowledge_full.tsv"
    def WHATIZ = "Tissue from https://tissues.jensenlab.org/Search  (${URL})"
"""
hostname 1>&2
mkdir -p TMP
export LC_ALL=C

jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP gtf2bed  --columns "gtf.feature,gene_name" -R "${fasta}"  "${gtf}" |\\
	awk -F '\t' '\$4=="gene" && \$5!="." && \$5!=""' |\\
	cut -f1,2,3,5 |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k4,4 -T TMP  |\\
	uniq > TMP/genes.bed

wget -O - "${URL}" |\\
    cut -f2,3,4 |\\
    awk -F '\t' '{OFS="\t";gsub(/[^A-Za-z0-9_:]+/,"_",\$3);print}' |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 |\\
    uniq > TMP/genes.txt

LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,1.2,1.3' TMP/genes.bed TMP/genes.txt |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq |\\
	bgzip > TMP/${TAG}.bed.gz

tabix --force -p bed TMP/${TAG}.bed.gz


mv TMP/*.bed.gz ./
mv TMP/*.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ}">' > ${TAG}.header
echo '##INFO=<ID=${TAG}_BRENDA,Number=.,Type=String,Description="${WHATIZ}">' >> ${TAG}.header


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
}



process ANNOTATE {
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(bed),path(bed_idx),path(header)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions

script:
    def TAG = task.ext.tag?:"TISSUES"
	def prefix=task.ext.prefix?:vcf.baseName+".tissues"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP OUTPUT



bcftools annotate \\
	--threads ${task.cpus} \\
	-a "${bed}" \\
	-h "${header}" \\
	--write-index \\
	-c "CHROM,POS,END,${TAG}_BRENDA,${TAG}" \\
	-O b \\
	--merge-logic '${TAG}:unique,${TAG}_BRENDA:unique' \\
	-o TMP/jeter.bcf \\
	'${vcf}'


mv TMP/jeter.bcf ${prefix}.bcf
mv TMP/jeter.bcf.csi ${prefix}.bcf.csi


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
