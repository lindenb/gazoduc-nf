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
		meta
		fasta
		fai
		dict
		bed
		vcfs /** tuple meta,vcf,vcf_index */
	main:

			annotate_ch = ANNOTATE(
				
				, vcfs)
			
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
tag "${meta.id?:vcf.name}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"

input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path("OTHER/*")
	tuple val(meta ),path(vcf),path(vcf_idx)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def ac= task.ext.max_ac?:10
	def TAG = task.ext.tag?:"ELSEWHERE"
	def prefix = task.ext.prefix?:vcf.baseName+".elsewhere"
"""
hostname 1>&2
mkdir -p TMP OUTPUT
set -x
find OTHER/ -name "*.vcf.gz" -o -name "*.bcf" |\\
	grep -v '\\.g\\.vcf\\.gz\$' |\\
	sort  > TMP/vcfs.list

echo -n '##INFO=<ID=${TAG},Number=.,Type=String,Description="Only for set for VCF having INFO/AC < ${ac}. Max displayed=${ac}. Samples found in other files: '>  TMP/${TAG}.header
cat TMP/vcfs.list| paste -sd ' ' | tr -d '\\n' >>   TMP/${TAG}.header
echo '">' >> TMP/${TAG}.header

# samples in that file
bcftools query -l "${vcf}" |\\
	sort -T TMP |\\
	uniq > TMP/samples.1.txt

# BED for this VCF
bcftools query -f "%CHROM\t%POS0\t%END\\n"  "${vcf}" |\\
	bedtools merge > TMP/roi.bed


#awk script for limit number of variants
cat << '_EOF_' > TMP/jeter.awk
BEGIN{
	OFS="\t";
	PREV="";
	N=0;
	}
	{
	if(NF>5) {
		G=\$NF;
		if(G=="0/0" || G=="./." || G=="0|0" || G==".|." || G=="0" || G==".") next;
		}

	S=sprintf("%s\t%s\t%s\t%s",\$1,\$2,\$3,\$4);
	if(S!=PREV) {
		if(N>0 && N<=${ac}) {
			for(i=1;i<=N;i++) {
				printf("%s\t%s\\n",PREV,samples[i]);
				}
			}
		N=0;
		PREV=S;
		}
	N++;
	if(N <= ${ac}) {
		samples[N] = \$5;
		}
	}
END {
	if(N>0 && N<=${ac}) {
		for(i=1;i<=N;i++) {
			printf("%s\t%s\\n",PREVS,samples[i]);
			}
		}
	}' 
_EOF_


cat  TMP/vcfs.list| while read F
do
	# samples in other file
	bcftools query -l "\${F}" |  sort -T TMP | uniq > TMP/samples.2.txt
	
	# adjust bed file
	jvarkit  -Djava.io.tmpdir=TMP  bedrenamechr -R "\${F}" --column 1 --convert SKIP  TMP/roi.bed > TMP/roi2.bed

	# remove common samples
	comm -13 TMP/samples.1.txt TMP/samples.2.txt > TMP/samples.3.txt

	# run query
	if test -s TMP/samples.3.txt && test -s TMP/roi2.bed
	then
		bcftools view --regions-file TMP/roi2.bed --samples-file TMP/samples.3.txt -O u "\${F}" |\\
		bcftools view -c 1 -O u |\\
		bcftools norm -f '${fasta}' --multiallelics -any -O u |\\
		bcftools query -i 'ALT!="*"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\\n]' |\\
		awk -F '\t' -f TMP/jeter.awk |\\
		jvarkit  -Djava.io.tmpdir=TMP  bedrenamechr -R "${fasta}" --column 1 --convert SKIP  |\\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4 > TMP/jeter2.bed

		# merge with previous bed if any
		if test -s TMP/jeter.bed
		then
			LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3 -k4,4  --merge TMP/jeter.bed TMP/jeter2.bed |\\
				uniq |\\
				awk -F '\t' -f TMP/jeter.awk  > TMP/jeter3.bed
			
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

	bcftools annotate \\
		--write-index \\
		--threads ${task.cpus} \\
		--keep-sites -i 'AC<${ac}' \\
		-a TMP/jeter.bed.gz \\
		--columns 'CHROM,POS,REF,ALT,${TAG}' \\
		--header-lines TMP/${TAG}.header \\
		--merge-logic '${TAG}:unique' \\
		-O b \\
		-o  TMP/${TAG}.bcf \\
		"${vcf}"
else

	bcftools view \\
		--write-index \\
		--threads ${task.cpus} \\
		-O b -o TMP/${prefix}.bcf "${vcf}"

fi

mv  TMP/${prefix}.bcf ./
mv  TMP/${prefix}.bcf.csi ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
}
