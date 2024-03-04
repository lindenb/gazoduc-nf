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
def TAG="CLINVAR"

workflow ANNOTATE_CLINVAR {
	take:
		genomeId
		vcfs /** json: vcf,vcf_index */
	main:

	
		if(hasFeature("clinvar") && !isBlank(params.genomes[genomeId],"clinvar_vcf_url")) {
			source_ch = DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.vcf, source_ch.index,vcfs)
			out1 = annotate_ch.output
			out2 = annotate_ch.count
			out3 = MAKE_DOC(genomeId).output
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
   	def url = genome.clinvar_vcf_url
    def whatis="ClinVar aggregates information about genomic variation and its relationship to human health."

"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p TMP

wget -O - "${url}" |\
	gunzip -c |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R "${genome.fasta}"  -n SKIP |\
	bcftools sort -T TMP/tmp -O b -o TMP/clinvar.bcf

bcftools view --header-only TMP/clinvar.bcf | grep "^##INFO" | cut -d '=' -f3 | cut -d ',' -f1 | grep -v '^CLN' | awk '{printf("INFO/%s\tCLN_%s\\n",\$1,\$1);}' > TMP/rename.tsv
bcftools annotate --rename-annots TMP/rename.tsv -O b -o TMP/jeter.bcf TMP/clinvar.bcf
mv TMP/jeter.bcf TMP/clinvar.bcf

bcftools index --force TMP/clinvar.bcf

mv TMP/clinvar.bcf ./${TAG}.db.bcf
mv TMP/clinvar.bcf.csi ./${TAG}.db.bcf.csi
"""
}


process MAKE_DOC {
executor "local"
input:
        val(genomeId)
output:
	path("${TAG}.html"),emit:output
script:
	def genome = params.genomes[genomeId]
   	def url = genome.clinvar_vcf_url
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>ClinVar aggregates information about genomic variation and its relationship to human health. <a href="${url}">${url}</a>.</dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	path(clinvar)
	path(clinvar_idx)
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


bcftools annotate -a "${clinvar}" -c "CLNSIG,CLN_ALLELEID" -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


###  
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > ${TAG}.count
mv TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
