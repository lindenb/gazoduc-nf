/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {getModules} from '../../modules/utils/functions.nf'

process VCF_TO_BED {
tag "${file(vcf).name}"
input:
	val(meta)
	val(vcf)
output:
	path("vcf2bed.bed"),emit:bed
	path("version.xml"),emit:version
script:

	if(vcf.endsWith(".list"))
	"""
	hostname 1>&2
	set -o pipefail
	module load ${getModules("bcftools")}

	cat "${vcf}" | while read V
	do
		bcftools index -s "\${V}" | awk -F '\t' -vV=\${V} '{printf("%s\t0\t%s\t%s\\n",\$1,\$2,V);}'
	done | sort -T . -t '\t' -k1,1 -k2,2n | uniq > vcf2bed.bed


	
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract contigs in multiple VCFs using bcftools</entry>
		<entry key="input">${vcf}</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF
	"""
	else
	"""
	hostname 1>&2
	set -o pipefail
	module load ${getModules("bcftools")}

	bcftools index -s "${vcf}" |\
		awk -F '\t' '{printf("%s\t0\t%s\t${vcf}\\n",\$1,\$2);}' |\
		sort -T . -t '\t' -k1,1 -k2,2n |\
		uniq > vcf2bed.bed


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract contigs in one VCFs using bcftools</entry>
		<entry key="input">${vcf}</entry>
		<entry key="bcftools">\$( bcftools --version-only)</entry>
	</properties>
	EOF
	"""
	}
