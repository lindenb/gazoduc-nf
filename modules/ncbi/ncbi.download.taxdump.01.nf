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
include {getVersionCmd;moduleLoad;runOnComplete} from '../utils/functions.nf'


process NCBI_DOWNLOAD_TAXDUMP_01 {
input:
	val(meta)
output:
	path("citations.dmp"),emit:citations
	path("delnodes.dmp"),emit:delnodes
	path("division.dmp"),emit:division
	path("gencode.dmp"),emit:gencode
	path("merged.dmp"),emit:merged
	path("names.dmp"),emit:names
	path("nodes.dmp"),emit:nodes
	path("gc.prt"),emit:gc
	path("version.xml"),emit:version
script:
	def url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
"""
hostname 1>&2
wget --no-check-certificate -O jeter.tar.gz "${url}"
tar xvfz jeter.tar.gz
rm jeter.tar.gz

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download ncbi taxdump</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}
