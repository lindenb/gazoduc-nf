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


def TAG="SPLICEAI"

 
workflow ANNOTATE_SPLICEAI {
	take:
		genomeId
		vcfs /** json vcf,vcf_index */
	main:
		if(hasFeature("spliceai") && !isBlank(params.genomes[genomeId],"gtf")) {
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
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>SpliceAI Prediction</dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
cpus 5
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
	def gtf = genome.gtf
	def db =genome.spliceai_annotation_type
	def distance = params.spliceai.splice_distance?:50
	def row = slurpJsonFile(json)	
"""
hostname 1>&2
${moduleLoad("bcftools bedtools htslib")}
mkdir -p TMP OUTPUT



tabix --regions "${row.bed}" "${gtf}" |\\
	awk -F '\t' '(\3=="exon") {for(i=0;i<2;i++) { {P=(i==0?int(\$4)-1:int(\$5));printf("%s\t%d\t%d\\n",\$1,P,P+1);}}' |\
	bedtools slop -b "${distance}" -g "${reference}.fai" |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n  > TMP/junctions.bed
	
if ! test -s TMP/junctions.bed
then
	tail -n1 "${reference}.fai"  '{printf("%s\t0\t1\\n",\$1);}' > TMP/junctions.bed
fi

bzip --force TMP/junctions.bed
tabix index --force -p bed TMP/junctions.bed.gz

bcftools view -O b --regions-file ^TMP/junctions.bed.gz -o TMP/off.bcf 
bcftools index --force TMP/off.bcf 

export OMP_NUM_THREADS=${task.cpus}

bcftools view --regions-file TMP/junctions.bed.gz |\\
	spliceai -R "${genome.fasta}"  -A "${db}" -D ${distance} |\
	bcftools view -O b -o TMP/in.bcf
	
bcftools index --force TMP/in.bcf 	

bcftools concat --allow-overlaps -O b -o  TMP/${TAG}.bcf  TMP/in.bcf  TMP/off.bcf
bcftools index --force TMP/${TAG}.bcf

cat << EOF > TMP/${TAG}.json
{
"vcf"   : "\${PWD}/OUTPUT/${TAG}.bcf",
"index" : "\${PWD}/OUTPUT/${TAG}.bcf.csi",
"bed"   : "${row.bed}"
}
EOF


bcftools query -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
