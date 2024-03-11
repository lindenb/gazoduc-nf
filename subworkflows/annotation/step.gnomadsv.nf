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
include {hasFeature;isBlank;backDelete;isHg19;isHg38} from './annot.functions.nf'


def TAG="GNOMADSV"

workflow ANNOTATE_GNOMADSV {
	take:
		genomeId
		vcfs /** json vcf,vcf_index */
	main:
		if(hasFeature("gnomadsv") && (isHg19(genomeId) || isHg38(genomeId)) ) {
			if(isHg19(genomeId)) {
				source_ch =  DOWNLOAD19(genomeId)
				}
			else if((isHg38(genomeId)) {
				source_ch =  DOWNLOAD38(genomeId)
				}
			else
				{
				throw new IllegalArgumentException("gnomad.sv:" + genomeId);
				}
			
			annotate_ch = ANNOTATE(source_ch.vcf, source_ch.index, vcfs)

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


process MAKE_DOC {
executor "local"
input:
        val(genomeId)
output:
	path("${TAG}.html"),emit:output
script:
	def genome = params.genomes[genomeId]
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>Frequent SV in Gnomad. min-AF:${params.annotations.gnomadsv_min_AF} Field ${params.annotations.gnomadsv_population} </dd>
</dl>
__EOF__
"""
}


process DOWNLOAD19 {
afterScript "rm -rf TMP"
memory "3G"
input:
	val(genomeId)
output:
	path("${TAG}.db.bcf"),emit:vcf
	path("${TAG}.db.bcf.csi"),emit:index
script:
	def genome = params.genomes[genomeId]
	def AF=params.annotations.gnomadsv_min_AF?:0.1
	def pop = params.annotations.gnomadsv_population?:"POPMAX_AF"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP


wget -O - "https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz" |\
	bcftools annotate --include 'INFO/${tag} >= ${af}' -x '^INFO/${tag}' -O v |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R "${genome.fasta}"  -n SKIP |\
	bcftools view -O b -o TMP/${TAG}.db.bcf

bcftools index --force TMP/${TAG}.db.bcf

mv TMP/${TAG}.* ./
"""
}


workflow DOWNLOAD38 {
take:
	genomeId
main:
	all_chr = Channel.of("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
	ch1 = DOWNLOAD38_CHR(genomeId,all_chr)
	ch2 = MERGE38(ch1.output.collect())
emit:
	vcf = ch2.vcf
	index = ch2.index
}


process DOWNLOAD38_CHR {
tag "${chr}"
afterScript "rm -rf TMP"
memory "3G"
input:
	val(genomeId)
	val(chr)
output:
	path("${TAG}.${chr}.db.bcf"),emit:output
script:
	def genome = params.genomes[genomeId]
	def AF=params.annotations.gnomadsv_min_AF?:0.1
	def pop = params.annotations.gnomadsv_population?:"POPMAX_AF"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
mkdir -p TMP

wget -O - "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/genome_sv/gnomad.v4.0.sv.chr${chr}.vcf.gz" |\
	bcftools annotate --include 'INFO/${tag} >= ${af}' -x '^INFO/${tag}' -O v |\
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfsetdict -R "${genome.fasta}"  -n SKIP |\
	bcftools view -O b -o TMP/${TAG}.${chr}.db.bcf

bcftools index --force TMP/${TAG}.${chr}.db.bcf

mv TMP/${TAG}.* ./
"""
}

process MERGE38 {
tag "${L}"
afterScript "rm -rf TMP"
input:
	val(L)
output:
	path("${TAG}.db.bcf"),emit:vcf
	path("${TAG}.db.bcf.csi"),emit:index
script:
	def pop = params.annotations.gnomadsv_population?:"POPMAX_AF"
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP


bcftools merge -a -O b  -o TMP/${TAG}.db.bcf ${L.join( " ")}
bcftools index --force TMP/${TAG}.db.bcf

mv TMP/${TAG}.* ./

"""
}


process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	path(vcf)
	path(index)
	path(json)
	//tuple path(vcf),path(vcf_idx),path(bed)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def pop = params.annotations.gnomadsv_population?:"POPMAX_AF"
	def AF=params.annotations.gnomadsv_min_AF?:0.1
	def whatis = "GNOMAD SV with ${pop} > ${AF}"
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools htslib")}
mkdir -p TMP OUTPUT

# check header exists
gunzip -c "${tabix}" | head -n 1 | tr "\t" "\\n" | grep -F '${pop}'

set -o pipefail

bcftools query --regions-file '${row.bed}' -f '%CHROM\t%POS0\t%END\t%ID|INFO/${pop}\\n' '${vcf}' > TMP/gnomad.bed
bgzip TMP/gnomad.bed
tabix --force -p bed TMP/gnomad.bed.gz

echo '##INFO=<ID=${TAG}_${pop},Number=.,Type=String,Description="Frequent SV in gnomad FORMAT: ID|AF . treshold AF>=:${AF}.">' >  TMP/${TAG}.header


bcftools annotate -a TMP/gnomad.bed.gz -c "CHROM,FROM,TO,${TAG}_${pop}" -h TMP/${TAG}.header  --merge-logic '${TAG}_${pop}:unique'  -O b -o TMP/${TAG}.bcf '${row.vcf}'
bcftools index TMP/${TAG}.bcf


cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT
${backDelete(row)}
"""
}
