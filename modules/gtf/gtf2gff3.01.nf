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

include {moduleLoad;isUrl;isBlank;getKeyValue;getModules;getBoolean} from '../utils/functions.nf'

process GTF_TO_GFF3_01 {
tag "${gtf.name}"
afterScript "rm -rf TMP"
memory "1g"
input:
	val(meta)
	path(gtf)
output:
	path("${gtf.getSimpleName()}.gff3${getBoolean(meta,"with_tabix")?".gz":""}"),emit:gff3
	path("${gtf.getSimpleName()}.gff3.gz.tbi"),emit:tbi,optional:true
	path("version.xml"),emit:version
script:
	def with_tabix = getBoolean(meta,"with_tabix")
	def concat = gtf.name.endsWith(".gz")?"gunzip -c":"cat"
	"""
	hostname 1>&2
	${moduleLoad("htslib gffread")}
	set -o pipefail
	mkdir TMP

	${concat} "${gtf}" |  gffread -FE |\
		LC_ALL=C sort -S ${task.memory.kilo} -T TMP -k1,1 -k4,4n > TMP/jeter.gff3

	if [ ! -z "${with_tabix?"Y":""}" ] ; then
		bgzip TMP/jeter.gff3
		tabix -p gff TMP/jeter.gff3.gz
		mv TMP/jeter.gff3.gz "${gtf.getSimpleName()}.gff3.gz"
		mv TMP/jeter.gff3.gz.tbi "${gtf.getSimpleName()}.gff3.gz.tbi"
	else
		mv TMP/jeter.gff3 "${gtf.getSimpleName()}.gff3"		
	fi
	
	#################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Convert GTF to GFF using gffread</entry>
		<entry key="gtf">${gtf}</entry>
		<entry key="with_tabix">${with_tabix}</entry>
		<entry key="gffread">\$(gffread --version 2>&1 )</entry>
		<entry key="bgzip">\$(bgzip --version | head -n 1)</entry>
		<entry key="tabix">\$(tabix --version | head -n 1)</entry>
	</properties>
	EOF
	"""
	
	stub:
	"""
	touch "${file(reference).getSimpleName()}.gff3${getBoolean(meta,"with_tabix")?".gz":""}"
        touch "${file(reference).getSimpleName()}.gff3.gz.tbi"
	echo "<properties/>" > version.xml
	"""
	}
