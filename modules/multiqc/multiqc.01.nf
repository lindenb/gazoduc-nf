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
include {moduleLoad;} from '../utils/functions.nf'

process MULTIQC_01 {
	tag "${row.files}"
	input:
		val(meta)
		val(row)
	output:
		path("${row.prefix?:""}multiqc.zip"),emit:zip
		path("${row.prefix?:""}multiqc/${row.prefix?:""}multiqc_report_data"),emit:datadir
		path("version.xml"),emit:version
	script:
		def prefix = row.prefix?:""
	"""
		hostname 1>&2
		${moduleLoad("multiqc")}

		mkdir -p "${prefix}multiqc"

		export LC_ALL=en_US.utf8
		multiqc  --filename  "${prefix}multiqc_report.html" \
			--title "${prefix}FASTQC"  \
			--force \
			--outdir "${prefix}multiqc" \
			--file-list "${row.files}"
		
		zip -9 -r "${prefix}multiqc.zip" "${prefix}multiqc"

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">MultiQC</entry>
	<entry key="wget.version">\$( multiqc --version )</entry>
</properties>
EOF


	"""

}


