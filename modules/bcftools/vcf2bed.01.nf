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
include {moduleLoad;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'


process VCF_TO_BED {
tag "${vcf.name}"
input:
	val(meta)
	path(vcf)
output:
	path("vcf2bed.bed"),emit:bed /* chrom start end vcf */
	path("contigs.txt"),emit:chromosomes /* uniq chromosome names */
	path("contigs.bed"),emit:chromsbed /* uniq chromosome names as BED */
	path("version.xml"),emit:version
script:
	if(meta.containsKey("with_header"))  throw new IllegalArgumentException("deprecated, use config/process");
	def with_header = task.ext.with_header?:false
	if(!(vcf.name.endsWith(".bed") || vcf.name.endsWith(".list") || vcf.name.endsWith(".vcf.gz") || vcf.name.endsWith(".bcf")))  throw new IllegalArgumentException("bad extension for vcf file:"+vcf);
	"""
	hostname 1>&2
	set -o pipefail
	${moduleLoad("bcftools")}

	if ${vcf.name.endsWith(".bed")} ; then

		grep -v '^#' '${vcf}' | cut -f1-4 | LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n | uniq > vcf2bed.bed

	elif ${vcf.name.endsWith(".list")} ; then
		cat "${vcf}" | while read V
		do
			bcftools index -s "\${V}" | awk -F '\t' -vV=\${V} '{printf("%s\t0\t%s\t%s\\n",\$1,\$2,V);}'
		done | LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n | uniq > vcf2bed.bed

	else

		bcftools index -s "${vcf.toRealPath()}" |\
			awk -F '\t' '{printf("%s\t0\t%s\t${vcf.toRealPath()}\\n",\$1,\$2);}' |\
			LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\
			uniq > vcf2bed.bed
	fi


	cut -f1 vcf2bed.bed | sort | uniq > contigs.txt 	
	cut -f1,2,3 vcf2bed.bed | sort | uniq > contigs.bed


	# add header
	if ${with_header} ; then

		echo -e 'contig\tstart\tend\tvcf' > jeter.tmp
		cat vcf2bed.bed >> jeter.tmp
		mv jeter.tmp vcf2bed.bed

		echo -e 'contig\tstart\tend' > jeter.tmp
		cat contigs.bed >> jeter.tmp
		mv jeter.tmp contigs.bed

	fi


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract contigs from VCF(s) using bcftools</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="version">${getVersionCmd("bcftools awk")}</entry>
	</properties>
	EOF
	"""
	}
