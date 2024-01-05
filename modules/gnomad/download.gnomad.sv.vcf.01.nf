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
include {isBlank;moduleLoad;isHg38;isHg19;getVersionCmd} from '../utils/functions.nf'

process DOWNLOAD_GNOMAD_SV_AS_VCF_01 {
input:
	val(meta)
	val(reference)
output:
	path("gnomad.sites${meta.suffix?:".bcf"}"),emit:vcf
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz":(isHg38(reference)?"":"")
	def suffix = meta.suffix?:".bcf"
	def tbi = suffix.contains("b")?"":"--tbi"
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail

if ! -z "${url}" ; then

wget -O - "${url}" |\
	bcftools view |\
	java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R "${reference}"  --onNotFound SKIP |\
	bcftools sort -T . -O ${suffix.contains("b")?"b":"z"} -o "gnomad.sites${suffix}"

else

cat << EOF > jeter.vcf
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
EOF

bcftools view -O ${suffix.contains("b")?"b":"z"} -o "gnomad.sites${suffix}" jeter.vcf
rm jeter.vcf

fi

bcftools index ${tbi} "gnomad.sites${suffix}"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download gnomad</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/vcfsetdict")}</entry>
        <entry key="versions"><a>${url}</a></entry>
</properties>
EOF
"""
}
