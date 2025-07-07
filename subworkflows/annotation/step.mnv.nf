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
include {moduleLoad;slurpJsonFile} from '../../modules/utils/functions.nf'
include {hasFeature;backDelete;isBlank} from './annot.functions.nf'
def TAG="MNV"

workflow ANNOTATE_MNV {
	take:
  		meta
		fasta
		fai
		dict
		gtf
		vcfs /** tuple vcf,vcf_index */
	main:

		if(hasFeature("mnv") && !isBlank(params.genomes[genomeId],"gtf")) {
			annotate_ch = ANNOTATE(genomeId, vcfs)
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
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>mnv with jvarkit</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
memory '3g'
input:
	val(genomeId)
	//tuple path(vcf),path(vcf_idx),path(bed)
	path(json)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
        def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools bedtools htslib jvarkit")}
mkdir -p TMP OUTPUT
set -x


bcftools norm -f '${reference}' --multiallelics -any -O u  '${row.vcf}'|\\
		bcftools view -i 'ALT!="*"' |\\
		java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfcombinetwosnvs -R '${reference}' -g '${genome.gtf}' > TMP/jeter1.vcf
		
bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/jeter1.bcf TMP/jeter1.vcf

rm  TMP/jeter1.bcf
bcftools index --force TMP/jeter1.bcf


bcftools annotate -a "TMP/jeter1.bcf" -c "CodonVariant" --merge-logic   'CodonVariant:unique'  -O b -o  TMP/${TAG}.bcf "${row.vcf}"
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
