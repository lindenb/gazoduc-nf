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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {hasFeature;isBlank;backDelete} from './annot.functions.nf'
def TAG="CONTRAST"

workflow ANNOTATE_CONTRAST {
	take:
		genomeId
		bed
		vcfs /** tuple vcf,vcf_index */
	main:

		if(hasFeature("contrast") && !isBlank(params,"pedigree") && !file(params.pedigree).name.equals("NO_FILE")) {
			annotate_ch = ANNOTATE(genomeId, file(params.pedigree),vcfs)
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
<dd>call <code>bcftools contrast</code> with <code>${params.pedigree}</code>.</dd>
</dl>
__EOF__
"""
}

process ANNOTATE {
tag "${json.name}"
afterScript "rm -rf TMP"
input:
	val(genomeId)
	path(pedigree)
	path(json)
	//tuple path(vcf),path(vcf_idx),path(bed)
output:
	//tuple path("OUTPUT/${TAG}.bcf"),path("OUTPUT/${TAG}.bcf.csi"),path(bed),emit:output
	path("OUTPUT/${TAG}.bcf"),emit:output
	path("OUTPUT/${TAG}.json"),emit:count
script:
	def args=""
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

# sample in VCF
bcftools query -l "${row.vcf}" | sort -T TMP | uniq > TMP/samples.txt

# samples in ped with phenotype
awk '{if(length(\$6) >0) print \$2}' '${pedigree}' | sort | uniq > TMP/jeter2.txt

# common samples
comm -12 TMP/samples.txt TMP/jeter2.txt > TMP/common.txt

# new pedigree
tr -s " " "\t" < '${pedigree}' | tr -s "\t" |\\
	sort -T TMP -t '\t' -k2,2 |\\
	join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6' - TMP/common.txt > TMP/jeter.ped


awk -F '\t' '{if(\$6=="case" || \$6=="2") print \$2;}' TMP/jeter.ped > TMP/cases.txt
awk -F '\t' '{if(\$6=="control" || \$6=="1") print \$2;}' TMP/jeter.ped > TMP/controls.txt

if test -s TMP/cases.txt && test -s TMP/controls.txt
then
	bcftools +contrast -0 TMP/controls.txt -1 TMP/cases.txt ${args} -O b -o  TMP/${TAG}.bcf 
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

##
bcftools query -N -f '.'  TMP/${TAG}.bcf | wc -c | awk '{printf("${TAG}\t%s\\n",\$1);}' > TMP/${TAG}.count

mv  TMP/${TAG}.* ./OUTPUT/
${backDelete(row)}
"""
}
