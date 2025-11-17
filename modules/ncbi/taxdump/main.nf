/*

Copyright (c) 2025 Pierre Lindenbaum

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


process NCBI_DOWNLOAD_TAXDUMP {
label "process_single"
input:
	val(meta)
output:
	tuple val(meta),path("citations.dmp"),emit:citations
	tuple val(meta),path("delnodes.dmp"),emit:delnodes
	tuple val(meta),path("division.dmp"),emit:division
	tuple val(meta),path("gencode.dmp"),emit:gencode
	tuple val(meta),path("merged.dmp"),emit:merged
	tuple val(meta),path("names.dmp"),emit:names
	tuple val(meta),path("nodes.dmp"),emit:nodes
	tuple val(meta),path("gc.prt"),emit:gc
	path("versions.xml"),emit:versions
script:
	def url = task.ext.url?:"https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
"""
curl -L -o jeter.tar.gz "${url}"
tar xvfz jeter.tar.gz
rm jeter.tar.gz

cat << EOF > versions.yml
"${task.process}"
    url: "${url}"
EOF
"""

stub:
"""
touch citations.dmp delnodes.dmp division.dmp gencode.dmp merged.dmp names.dmp nodes.dmp gc.prt versions.yml
"""
}
