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


def TAG="GENCC"
def URL = "http://ftp.ebi.ac.uk/pub/databases/GO/goa/bhf-ucl/gene_association.goa_bhf-ucl.gz"
def WHATIZ = "gencc The Gene Curation Coalition brings together groups engaged in the evaluation of gene-disease validity with a willingness to share data publicly https://thegencc.org/about.html"
 
workflow ANNOTATE_GENCC {
	take:
		genomeId
		vcfs /** json vcf,vcf_index */
	main:
		if(hasFeature("gencc") && !isBlank(params.genomes[genomeId],"gtf")) {
			source_ch =  DOWNLOAD(genomeId)
			annotate_ch = ANNOTATE(source_ch.bed, source_ch.tbi,source_ch.header,vcfs)

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

process DOWNLOAD{
afterScript "rm -rf TMP"
memory "2g"
input:
	val(genomeId)
output:
	path("${TAG}.bed.gz"),emit:bed
	path("${TAG}.bed.gz.tbi"),emit:tbi
	path("${TAG}.header"),emit:header
script:
	def genome = params.genomes[genomeId]
	
    def reference = genome.fasta
"""

hostname 1>&2
${moduleLoad("htslib bedtools jvarkit")}
set -o pipefail
mkdir -p TMP


java -jar \${JVARKIT_DIST}/jvarkit.jar gtf2bed --columns "gtf.feature,gene_name" -R "${reference}" "${genome.gtf}" |\\
	awk -F '\t' '(\$4=="gene")' |\\
	cut -f1,2,3,5 |\\
	LC_ALL=C sort -T TMP -t '\t' -k4,4 |\\
	uniq > TMP/genes.bed

wget --no-check-certificate -q  -O - "${URL}" |\\
	awk -F '\t' '(NR==1) {C=-1;D=-1;for(i=1;i<=NF;i++) {if(\$i=="\\"gene_symbol\\"") C=i;if(\$i=="\\"disease_title\\"") D=i;}next;} {if(C<1 ||  D<1) next; G=\$C;H=\$D; gsub(/[^A-Za-z0-9\\.\\-]+/,"_",G);gsub(/[^A-Za-z0-9\\.\\-]+/,"_",H);  printf("%s\t%s\\n",G,H);}'  |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1  | uniq > TMP/org.tsv


LC_ALL=C  join -t '\t' -1 4 -2 1 -o '1.1,1.2,1.3,2.2' TMP/genes.bed TMP/org.tsv |\\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bgzip > TMP/${TAG}.bed.gz


tabix --force -p bed TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${WHATIZ} ${URL}.">' > ${TAG}.header
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
	def url = genome.greendb_url
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>${WHATIZ} <a href="${URL}">${URL}</a></dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	path(tabix)
	path(tbi)
	path(header)
	//tuple path(vcf),path(vcf_idx),path(bed)
	path(json)
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

bcftools annotate -a "${tabix}" -h "${header}" -c "CHROM,FROM,TO,${TAG}"  --merge-logic '${TAG}:unique'  -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
