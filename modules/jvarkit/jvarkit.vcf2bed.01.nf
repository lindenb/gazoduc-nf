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

include {moduleLoad;getBoolean;assertNotEmpty;getKeyValue} from '../utils/functions.nf'

process VCF_TO_BED_01 {
	tag "${file(vcf).name}"
	memory "2g"
	input:
		val(meta)
		val(vcf)
	output:
		path("bed2vcf.tsv"),emit:output
		path("version.xml"),emit:version
	script:
		def distance = meta.vcf2interval_distance?:"10mb"
		def min_distance = meta.vcf2interval_min_distance?:"1kb"
	"""
	hostname 1>&2
	${moduleLoad("jvarkit bcftools")}
	set -o pipefail

	mkdir BEDS

	bcftools view -G "${vcf}" |\
	java -jar \${JVARKIT_DIST}/vcf2intervals.jar \
		--bed \
		--distance "${distance}" \
		--min-distance "${min_distance}" |\
	cut -f1,2,3 |\
	split -a 9 --lines=1 --additional-suffix=.bed - "BEDS/cluster."

	find \${PWD}/BEDS -type f -name "*.bed" |\
		awk  'BEGIN{printf("bed\tvcf\\n");}{printf("%s\t${vcf}\\n",\$1);}' > bed2vcf.tsv
	
	#######################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">convert vcf to bed using jvarkit/vcf2intervals</entry>
		<entry key="vcf">${vcf}</entry>
		<entry key="distance">${distance}</entry>
		<entry key="min_distance">${min_distance}</entry>
		<entry key="vcf2interval.version">\$(java -jar \${JVARKIT_DIST}/vcf2intervals.jar  --version)</entry>
	</properties>
	EOF
	"""
	stub:
	"""
	touch bed2vcf.tsv
	echo "<properties/>" > version.xml
	"""
	}
