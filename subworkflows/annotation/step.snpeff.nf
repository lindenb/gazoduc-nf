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

def TAG="SNPEFF"

workflow ANNOTATE_SNPEFF {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:
		if(hasFeature("snpeff") && !isBlank(params.genomes[genomeId],"gtf")) {
			source_ch =  BUILD_SNPEFF(genomeId)
			annotate_ch = ANNOTATE(genomeId,source_ch.output, vcfs)
			
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

/** build snpeff Database from gtf */
process BUILD_SNPEFF {

afterScript "rm -f org.tsv genes.tsv snpeff.errors"
memory "10g"
input:
        val(genomeId)
output:
       	path("snpEff.config"),emit:output
script:
        def genome = params.genomes[genomeId]
        def dbName = genomeId
        def reference = genome.fasta
        def gtf = genome.gtf
"""
hostname 1>&2
${moduleLoad("snpEff")}
set -o pipefail

mkdir -p "data/${dbName}"
ln -s "${reference}" "data/${dbName}/sequences.fa"

cp "${gtf}"  data/${dbName}/genes.gtf.gz
gunzip data/${dbName}/genes.gtf.gz

# write snpEff contig
cat << EOF > snpEff.config
data.dir = \${PWD}/data/
${dbName}.genome = Human
${dbName}.reference = ${reference}
EOF

# build database
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${SNPEFF_JAR}  build -gtf22 -v "${dbName}" 2> snpeff.errors

rm snpeff.errors

test -s "data/${dbName}/snpEffectPredictor.bin"

rm  data/${dbName}/genes.gtf
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
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>Annotation with SNPEFF : ${genome.gtf}</dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
memory '3g'
input:
	val(genomeId)
	path(config)
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
${moduleLoad("snpEff bcftools")}
mkdir -p TMP OUTPUT

set -o pipefail

bcftools view '${row.vcf}' -O v |\
java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar "\${SNPEFF_JAR}" eff -config '${config}' \
                                -nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf '${genomeId}' > TMP/jeter1.vcf

bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/${TAG}.bcf TMP/jeter1.vcf
bcftools index TMP/${TAG}.bcf
rm TMP/jeter1.vcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


###
bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
