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
include {hasFeature;backDelete} from './annot.functions.nf'
def TAG="ELSEWHERE"

workflow ANNOTATE_ELSEWHERE {
	take:
		genomeId
		bed
		vcfs /** tuple vcf,vcf_index */
	main:

		if(hasFeature("elsewhere") && !params.annotations.elsewhere.list.equals("NO_FILE")) {
			annotate_ch = ANNOTATE(genomeId, vcfs)
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

process MAKE_DOC {
executor "local"
output:
	path("${TAG}.html"),emit:output
script:
"""
cat << __EOF__ > ${TAG}.html
<dl>
<dt>${TAG}</dt>
<dd>Sample with ALT genotypes found in the following file(s):<ul>
__EOF__

grep -v "#"  '${params.annotations.elsewhere.list}' | awk '{printf("<li>%s</li>\\n",\$0);}'  >> ${TAG}.html

cat << __EOF__ >> ${TAG}.html
</ul>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
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
	def others = params.annotations.elsewhere.list
	def ac= (params.annotations.elsewhere.max_ac?:10)
"""
hostname 1>&2
${moduleLoad("bcftools bedtools htslib")}
mkdir -p TMP OUTPUT
set -x

echo -n '##INFO=<ID=${TAG},Number=.,Type=String,Description="Only for set for VCF having INFO/AC < ${ac}. Max displayed=${ac}. Samples found in other files: '>  TMP/${TAG}.header
sort "${others}" | paste -sd ' ' | tr -d '\\n' >>   TMP/${TAG}.header
echo '">' >> TMP/${TAG}.header


bcftools query -l "${row.vcf}" | sort -T TMP | uniq > TMP/samples.1.txt
bcftools query -f "%CHROM\t%POS0\t%END\\n"  "${row.vcf}" | bedtools merge > TMP/roi.bed

grep -v '^#' "${others}" | grep -E '\\.(vcf\\.gz|bcf)\$' | grep -v '\\.g\\.vcf\\.gz\$' | while read F
do
	# samples in other file
	bcftools query -l "\${F}" |  sort -T TMP | uniq > TMP/samples.2.txt
	
	# remove common samples
	comm -13 TMP/samples.1.txt TMP/samples.2.txt > TMP/samples.3.txt
	if test -s TMP/samples.3.txt && test -s TMP/roi.bed
	then
		bcftools view --regions-file TMP/roi.bed --samples-file TMP/samples.3.txt -O u "\${F}" |\\
		bcftools view -c 1 -O u |\\
		bcftools norm -f '${reference}' --multiallelics -any -O u |\\
		bcftools query -i 'ALT!="*"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\\n]' |\\
		awk -F '\t' '{G=\$NF;if(G ~ /^[0\\.][|/][0\\.]\$/ || G=="0" || G==".") next; print;}' |\\
		awk -F '\t' 'BEGIN{PREV="";N=0;} {S=sprintf("%s:%s:%s:%s",\$1,\$2,\$3,\$4);if(S!=PREV) {print;N=0;PREV=S;} else if(N< ${ac}) {print;} N++;}' |\\
		cut -f 1-5 |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4 > TMP/jeter2.bed

		# merge with previous bed if any
		if test -s TMP/jeter.bed
		then
			LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4  --merge TMP/jeter.bed TMP/jeter2.bed | uniq |\\
			awk -F '\t' 'BEGIN{PREV="";N=0;} {S=sprintf("%s:%s:%s:%s",\$1,\$2,\$3,\$4);if(S!=PREV) {print;N=0;PREV=S;} else if(N< ${ac}) {print;} N++;}' > TMP/jeter3.bed
			mv TMP/jeter3.bed TMP/jeter.bed
			rm TMP/jeter2.bed
		else
			mv TMP/jeter2.bed TMP/jeter.bed
		fi	
	fi
done

if test -s TMP/jeter.bed
then
	bgzip TMP/jeter.bed
	tabix --force -s 1 -b 2 -e 2 TMP/jeter.bed.gz
	bcftools annotate --keep-sites -i 'AC<${ac}' -a TMP/jeter.bed.gz --columns 'CHROM,POS,REF,ALT,${TAG}' --header-lines TMP/${TAG}.header --merge-logic '${TAG}:unique' -O b -o  TMP/${TAG}.bcf "${row.vcf}"
else

	bcftools view -O b -o TMP/${TAG}.bcf "${row.vcf}"

fi

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
