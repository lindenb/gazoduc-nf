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
include { paramsToString; escapeXml } from './functions.nf'

process PARAMS_MULTIQC {
executor "local"
input:
      	val(meta)
output:
        path("params_mqc.html"),emit:output
	path("version.xml"),emit:version
script:
	def str = escapeXml(paramsToString(params));
        """

cat << '__EOF__' > params_mqc.html
<!--
id: 'nf_params'
section_name: 'Nextflow Paramaters'
description: 'Parameters used in nextflow.'
-->
<pre>${str}</pre>
__EOF__


	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="Name">${task.process}</entry>
		<entry key="Description">Dump params as MULTIQC</entry>
	</properties>
	EOF
        """
	}

