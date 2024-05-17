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

include {slurpJsonFile;parseBoolean;moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete;isHg38} from './annot.functions.nf'

def TAG="UNIPROT"

workflow ANNOTATE_UNIPROT {
	take:
		genomeId
		bed
		vcfs /** json: vcf,vcf_index */
	main:

             	if(hasFeature("uniprot") && isHg38(genomeId) ) {
                        features_ch = Channel.of("zn_fing","binding","coiled","disulfide","metal","motif")
                        source_ch = DOWNLOAD(genomeId,features_ch)
                        annotate_ch = ANNOTATE(source_ch.output.map{T->T.join("\t")}.collect(),vcfs)
                        out1 = annotate_ch.output
                        out2 = annotate_ch.count
			            out3 =  MAKE_DOC(genomeId).output
                        }
                else
                    	{
                        out1 = vcfs
                        out2 = Channel.empty()
                        out3 = Channel.empty()
                        }

	emit:
		output = annotate_ch.output
		count = annotate_ch.count
		doc = out3
}

process DOWNLOAD{
tag "${type}"
afterScript "rm -rf TMP"
memory "2g"
input:
	val(genomeId)
    val(type)
output:
	tuple val("${TAG}_${type}"),path("${TAG}.${type}.bed.gz"),path("${TAG}.${type}.bed.gz.tbi"),path("${TAG}.${type}.header"),emit:output
script:
	def genome = params.genomes[genomeId]
	def url = "https://trackhubs.uniprot.org/genome_annotation_tracks/UP000005640_9606_hub/hg38/UP000005640_9606_${type}.bb"
    def reference = genome.fasta

"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("htslib jvarkit bedtools ucsc")}

set -o pipefail
wget -O TMP/jeter.bb "${url}"

bigBedToBed TMP/jeter.bb stdout |\\
        cut -f1-3 |\\
        java -jar \${JVARKIT_DIST}/jvarkit.jar bedrenamechr -f "${reference}" --column 1 --convert SKIP  |\\
        LC_ALL=C sort -t '\t' -k1,1 -k2,2n -T TMP |\\
        bedtools merge |\\
        sed 's/\$/\t1/' > TMP/${TAG}.${type}.bed.gz

tabix -p bed -f TMP/${TAG}.${type}.bed.gz

mv TMP/${TAG}.${type}.bed.gz ./
mv TMP/${TAG}.${type}.bed.gz.tbi ./

echo '##INFO=<ID=${TAG}_${type},Number=0,Type=Flag,Description="feature ${type} in uniprot ">' > ${TAG}.${type}.header
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
	def url = genome.vista_enhancers_url
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>VISTA enhancers. <a href="${url}">${url}</a></dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	val(L)
	path(json)
output:
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def  row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

bcftools view -b -O b -o TMP/jeter.bcf ${row.vcf}

cat TMP/jeter.list while read COL BED TBI HEADER
do
    bcftools annotate -a "\${BED}" -h "\${HEADER}" -c "CHROM,FROM,TO,\${COL}" -O b -o TMP/jeter2.bcf  TMP/jeter.bcf
    mv TMP/jeter2.bcf TMP/jeter.bcf
done

mv TMP/jeter.bcf  TMP/${TAG}.bcf
bcftools index --force TMP/${TAG}.bcf



cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
