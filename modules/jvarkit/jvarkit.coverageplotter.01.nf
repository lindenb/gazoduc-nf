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

process JVARKIT_COVERAGE_PLOTTER_01 {
	tag "${row.interval}"
	memory "5g"
	input:
		val(meta)
		val(row)
	output:
		tuple row,path("clusters.list"),emit:output
		path("version.xml"),emit:version
	script:
		def interval = row.interval
		def reference = row.reference
		def gtf = row.gtf?:""
		def known = row.known?:""
		def output = (meta.prefix?:"")+interval.replaceAll("[:-]+","_")+".html"
	"""
	hostname 1>&2
	set -o pipefail
	${moduleLoad("jvarkit")}

	test ! -z "${row.interval}"

	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${JVARKIT_DIST}/coverageplotter.jar \
		-R "${reference}" \
		${isBlank(gtf)?"":"--gtf \"${gtf}\"} \
		${isBlank(known)?"":"--known \"${known}\"} \
		--region "${interval}" \
		${row.bams} > ${output}

	#####################################################################################

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">plot normalized coverage for bams</entry>
		<entry key="interval">${interval}</entry>
		<entry key="reference">${reference}</entry>
		<entry key="gtf">${gtf}</entry>
		<entry key="method">${method}</entry>
		<entry key="coverageplotter.version">\$(java  -jar \${JVARKIT_DIST}/coverageplotter.jar --version)</entry>
	</properties>
	EOF
	"""
	stub:
	"""
	touch "${output}"
	echo "<properties/>" > version.xml
	"""
	}
