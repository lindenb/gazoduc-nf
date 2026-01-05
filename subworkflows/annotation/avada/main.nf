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
include {isGRCh37} from '../../../modules/utils/k1.nf'


workflow AVADA {
	take:
		meta
		fasta
		fai
		dict
		vcfs /* meta, vcf,vcf_index */
	main:
	
		if(isGRCh37(fai) ) {
			DOWNLOAD(fasta,fai,dict)
			ANNOTATE(DOWNLOAD.out.vcf,vcfs)
			vcf = ANNOTATE.out.vcf
			doc = ANNOTATE.out.doc
			versions = ANNOTATE.out.versions
			}
		else
			{
			doc = Channel.empty()
			versions = Channel.empty()
			vcf = vcfs
			}
	emit:
		doc
		versions
		vcf
	}


process DOWNLOAD {
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	path val(meta1),path("*.db.bcf"),path("*.db.bcf.csi"),emit:vcf
	path("versions.yml"),emit:versions
	path("doc.md"),emit:doc
script:
	def TAG=task.ext.tag?:"AVADA"
	def url = task.ext.url?:"http://bejerano.stanford.edu/AVADA/avada_v1.00_2016.vcf.gz"
"""
hostname 1>&2
mkdir -p TMP

wget -O - "${url}" |\\
	gunzip -c |\\
	awk -F '\t' '/^#CHROM/ {printf("##INFO=<ID=PMID,Number=.,Type=String,Description=\\"AVADA PMID\\">\\n##INFO=<ID=GENE_SYMBOL,Number=.,Type=String,Description=\\"AVADA GENE_SYMBOL\\">\\n##INFO=<ID=ENSEMBL_ID,Number=.,Type=String,Description=\\"AVADA ENSEMBL_ID\\">\\n##INFO=<ID=ENTREZ_ID,Number=.,Type=String,Description=\\"AVADA ENTREZ_ID\\">\\n##INFO=<ID=REFSEQ_ID,Number=.,Type=String,Description=\\"AVADA REFSEQ_ID\\">\\n##INFO=<ID=STRAND,Number=.,Type=String,Description=\\"AVADA STRAND\\">\\n##INFO=<ID=ORIGINAL_VARIANT_STRING,Number=.,Type=String,Description=\\"AVADA ORIGINAL_VARIANT_STRING\\">\\n\\n");} {print;}' |\\
	sed 's/PMID/${TAG}_PMID/g' |\
	jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP vcfsetdict -R "${fasta}"  -n SKIP |\\
	bcftools sort -T TMP/tmp -O b -o TMP/${TAG}.db.bcf

bcftools view --header-only TMP/${TAG}.db.bcf | grep "^##INFO" | cut -d '=' -f3 | cut -d ',' -f1 | grep -v '^PMID' | awk '{printf("INFO/%s\t${TAG}_%s\\n",\$1,\$1);}' > TMP/rename.tsv
bcftools annotate --rename-annots TMP/rename.tsv -O b -o TMP/jeter.bcf TMP/${TAG}.db.bcf
mv TMP/jeter.bcf TMP/${TAG}.db.bcf

bcftools index --threads '${task.cpus}' --force TMP/${TAG}.db.bcf

mv TMP/${TAG}.db.bcf ./
mv TMP/${TAG}.db.bcf.csi ./



cat << 'EOF' > doc.md
# annotations:avada

`INFO/${TAG}`

> pubmed-id in avada. The AVADA database includes **unvalidated** variant evidence data, 
> automatically retrieved from 61,116 full text papers deposited in PubMed until 07-2016

EOF

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}


process ANNOTATE {
tag "${meta.id?vcf.name}"
afterScript "rm -rf TMP"
input:
	tuple val(met1),path(database),path(database_idx)
	tuple val(meta),path(vcf),path(vcf_idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def TAG=task.ext.tag?:"AVADA"
"""
hostname 1>&2
mkdir -p TMP OUTPUT

bcftools annotate \\
	--write-index \\
	--threads ${task.cpus} \\
	-a "${database}" \\
	-c "${TAG}_PMID" \\
	--merge-logic '${TAG}_PMID:unique' \\
	-O b \\
	-o TMP/${TAG}.bcf \\
	'${vcf}'


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}