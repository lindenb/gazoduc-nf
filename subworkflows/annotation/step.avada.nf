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
include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete} from './annot.functions.nf'
def TAG="AVADA"

workflow ANNOTATE_AVADA {
	take:
		genomeId
		vcfs /** json : tuple vcf,vcf_index */
	main:

	
		if(hasFeature("avada") && !isBlank(params.genomes[genomeId],"ucsc_name") && params.genomes[genomeId].ucsc_name.equals("hg19") ) {
			source_ch = DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.vcf, source_ch.index,vcfs)
			out1 = annotate_ch.output
			out2 = annotate_ch.count
			out3 = MAKE_DOC().output
			}
		else
			{
			out1 = vcfs
			out2 = Channel.empty()
			out3 = Channel.empty()
			}
	emit:
		output = out1
		count = out2
		doc = out3
	}


process DOWNLOAD {
afterScript "rm -rf TMP"
memory "3g"
input:
        val(genomeId)
output:
       	path("${TAG}.db.bcf"),emit:vcf
        path("${TAG}.db.bcf.csi"),emit:index
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
       	def url = "http://bejerano.stanford.edu/AVADA/avada_v1.00_2016.vcf.gz"
	def whatis = "pubmed-id in avada. The AVADA database includes unvalidated variant evidence data, automatically retrieved from 61,116 full text papers deposited in PubMed until 07-2016"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\\
	gunzip -c |\\
	awk -F '\t' '/^#CHROM/ {printf("##INFO=<ID=PMID,Number=.,Type=String,Description=\\"AVADA PMID\\">\\n##INFO=<ID=GENE_SYMBOL,Number=.,Type=String,Description=\\"AVADA GENE_SYMBOL\\">\\n##INFO=<ID=ENSEMBL_ID,Number=.,Type=String,Description=\\"AVADA ENSEMBL_ID\\">\\n##INFO=<ID=ENTREZ_ID,Number=.,Type=String,Description=\\"AVADA ENTREZ_ID\\">\\n##INFO=<ID=REFSEQ_ID,Number=.,Type=String,Description=\\"AVADA REFSEQ_ID\\">\\n##INFO=<ID=STRAND,Number=.,Type=String,Description=\\"AVADA STRAND\\">\\n##INFO=<ID=ORIGINAL_VARIANT_STRING,Number=.,Type=String,Description=\\"AVADA ORIGINAL_VARIANT_STRING\\">\\n\\n");} {print;}' |\\
	sed 's/PMID/${TAG}_PMID/g' |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R "${reference}"  -n SKIP |\\
	bcftools sort -T TMP/tmp -O b -o TMP/${TAG}.db.bcf

bcftools view --header-only TMP/${TAG}.db.bcf | grep "^##INFO" | cut -d '=' -f3 | cut -d ',' -f1 | grep -v '^PMID' | awk '{printf("INFO/%s\t${TAG}_%s\\n",\$1,\$1);}' > TMP/rename.tsv
bcftools annotate --rename-annots TMP/rename.tsv -O b -o TMP/jeter.bcf TMP/${TAG}.db.bcf
mv TMP/jeter.bcf TMP/${TAG}.db.bcf

bcftools index --force TMP/${TAG}.db.bcf

mv TMP/${TAG}.db.bcf ./
mv TMP/${TAG}.db.bcf.csi ./
"""
}


process MAKE_DOC {
executor "local"
output:
	path("${TAG}.html"),emit:output
script:
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>pubmed-id in avada. The AVADA database includes <b>unvalidated</b> variant evidence data, automatically retrieved from 61,116 full text papers deposited in PubMed until 07-2016</dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	path(database)
	path(database_idx)
	path(json)
	//tuple path(vcf),path(vcf_idx),path(bed)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	 def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

bcftools annotate -a "${database}" -c "${TAG}_PMID" --merge-logic '${TAG}_PMID:unique' -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index --force TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
