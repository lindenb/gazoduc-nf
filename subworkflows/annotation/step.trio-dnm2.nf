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


def TAG="DNM2"

workflow ANNOTATE_DNM2 {
	take:
		genomeId
		vcfs /** tuple vcf,vcf_index */
	main:

		if(hasFeature("dbnm2") && !isBlank(params.pedigree)) {
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
<dd>call <code>bcftools  +trio-dnm2</code> with <code>${params.pedigree}</code>.</dd>
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
	path("OUTPUT/${TAG}.json"),emit:output
	path("OUTPUT/${TAG}.count"),emit:count
script:
	def args="--use-NAIVE "
	def row = slurpJsonFile(json)
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUTPUT

# sample in VCF
# new pedigree
tr -s " " "\t" < '${pedigree}' | tr -s "\t" |\\
	cut -f1-4 |\
	sort -T TMP -t '\t' -k2,2 > TMP/tmp.ped

# get samples with parent
awk -F '\t' '(\$3!="0" && \$4!="0")' TMP/tmp.ped | cut -f2,3,4 | tr "\t" "\n" | sort -T TMP | uniq > TMP/jeter.txt

# common samples
join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4' TMP/tmp.ped TMP/jeter.txt  > TMP/tmp2.ped
mv TMP/tmp2.ped TMP/tmp.ped


if test -s TMP/tmp.ped
then
	bcftools +trio-dnm2 -P TMP/tmp.ped ${args} -O b -o  TMP/${TAG}.bcf 
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
mv -v TMP/${TAG}.* OUTPUT/
${backDelete(row)}
"""
}
