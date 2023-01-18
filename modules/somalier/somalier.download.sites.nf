/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {moduleLoad;assertKnownReference;getKeyValue;isHg19;isHg38} from '../utils/functions.nf'

process SOMALIER_DOWNLOAD_SITES {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("sites.vcf.gz"), emit: vcf
	path("sites.vcf.gz.tbi"), emit:index
	path("version.xml"), emit:version
script:
	def fasta = assertKnownReference(reference)
	def url=(isHg19(fasta)?"https://github.com/brentp/somalier/files/3412453/sites.hg19.vcf.gz":
		(isHg38(fasta)?"https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz":
		""))
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail
set -x

wget -O - "${url}" |\
	gunzip -c |\
	java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R "${reference}"  --onNotFound SKIP |\
	bcftools sort -T . -o sites.vcf.gz -O z

bcftools index -t sites.vcf.gz

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download VCF sites for somalier</entry>
        <entry key="url"><a>${url}</a></entry>
	<entry key="bcftools.version">\$( bcftools --version-only)</entry>
</properties>
EOF
"""
stub:
"""
touch sites.vcf.gz sites.vcf.gz.tbi
echo "<properties/>" > version.xml
"""
}
