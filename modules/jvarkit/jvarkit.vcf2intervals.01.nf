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

include {moduleLoad;getVersionCmd;parseBoolean} from '../utils/functions.nf'

process JVARKIT_VCF_TO_INTERVALS_01 {
	tag "${row.vcf.name} ${row.contig?:""} ${row.interval?:""}"
	memory "2g"
	input:
		val(meta)
		val(row) /* map : vcf , contig */
	output:
		path("intervals.bed"),emit:bed /* chrom/start/end/vcf */
		path("version.xml"),emit:version
	script:
		def contig = row.contig?:""
		def interval = row.interval?:""
		def vcf = row.vcf
		def with_header = parseBoolean(row.with_header)
		def distance = meta.vcf2interval_distance?:"10mb"
		def min_distance = meta.vcf2interval_min_distance?:"1kb"
	"""
	hostname 1>&2
	${moduleLoad("jvarkit bcftools")}
	set -o pipefail

	mkdir BEDS

	if ${with_header} ; then
		echo -e 'contig\tstart\tend\tvcf' > intervals.bed
	fi


	bcftools view -G "${vcf.toRealPath()}" ${interval.isEmpty()?"":"\"${interval}\""}  ${contig.isEmpty()?"":"\"${contig}\""} |\
	java -jar \${JVARKIT_DIST}/vcf2intervals.jar \
		--bed \
		--distance "${distance}" \
		--min-distance "${min_distance}" |\
		awk  '{printf("%s\t%s\t%s\t${vcf.toRealPath()}\\n",\$1,\$2,\$3);}' >> intervals.bed
	
	#######################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">convert vcf to bed using jvarkit/vcf2intervals</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="contig">${contig}</entry>
		<entry key="interval">${interval}</entry>
		<entry key="distance">${distance}</entry>
		<entry key="min_distance">${min_distance}</entry>
		<entry key="versions">${getVersionCmd("jvarkit/vcf2intervals awk bcftools")}</entry>
	</properties>
	EOF
	"""
	}
