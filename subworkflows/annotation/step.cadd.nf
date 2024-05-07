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

def TAG="CADD"

String getUrl(genomeId) {
        if(isHg19(genomeId)) return "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh37/whole_genome_SNVs.tsv.gz";
	if(isHg38(genomeId)) return "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
        return "";
        }


workflow ANNOTATE_CADD {
	take:
		genomeId
		vcfs /** json: vcf,index,bed */
	main:

		 if(hasFeature("cadd") && (isHg19(genomeId) || isHg38(genomeId)) ) {
			    annotate_ch = ANNOTATE(genomeId,vcfs)
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
	def url = getUrl(genomeId)
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>CADD annotation from <a href="${url}">${url}</a></dd>
</dl>
__EOF__
"""
}


process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP *.tsv.gztbi"
memory "5g"
maxForks 1
input:
	val(genomeId)
	path(json)
output:
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def url = getUrl(genomeId)
	def genome = params.genomes[genomeId]
	def row = slurpJsonFile(json)

"""
hostname 1>&2
${moduleLoad("htslib bcftools jvarkit bedtools")}
mkdir -p TMP OUTPUT
set -o pipefail

# get intervals for this VCF, remove chr prefix, sort and merge
bcftools query -f '%CHROM\t%POS0\t%END\\n' '${row.vcf}' |\\
    sed 's/^chr//' |\\
    LC_ALL=C sort -T TMP -t $'\t' -k1,1 -k2,2n |\\
    bedtools merge > TMP/intervals.bed

tabix -h  \\
    --regions TMP/intervals.bed \\
    "${url}" |\\
    java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${genome.fasta}" --column 1 --convert SKIP  |\\
    bgzip > TMP/database.tsv.gz

tabix -c 1 -c '#' -b 2 -e 2 TMP/database.tsv.gz


bcftools annotate -c 'CHROM,POS,REF,ALT,${TAG}_RAWSCORE,${TAG}_PHRED' \\
    -a TMP/database.tsv.gz \\
    -H '##INFO=<ID=${TAG}_RAWSCORE,Number=1,Type=Float,Description="Raw Score in CADD ${url}">' \\
    -H '##INFO=<ID=${TAG}_PHRED,Number=1,Type=Float,Description="Phred Score in CADD ${url}">' \\
    -O b -o  TMP/${TAG}.bcf \\	
    '${row.vcf}'


bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF

###
bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
