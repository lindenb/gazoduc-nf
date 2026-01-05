/*

Copyright (c) 2026 Pierre Lindenbaum

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

include {escapeXml} from './functions.nf'

params.downstream_cmd = " false"

process CONCAT_FILES_02 {
tag "N=${L.size()}"
executor "local"
input:
        val(L)
output:
        path("concat${params.suffix?:".list"}"),emit:output
	path("version.xml"),emit:version
script:
	def concat_n_files = params.concat_n_files?:"50"
	def downstream_cmd = params.downstream_cmd?:""
"""
hostname 1>&2
cat << EOF | xargs -L ${concat_n_files} cat ${downstream_cmd} > "concat${params.suffix?:".list"}"
${L.join("\n")}
EOF

sleep 5
#############################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="Name">${task.process}</entry>
        <entry key="Description">Concatenate the content of N='${L.size()}' file(s)</entry>
        <entry key="downstream.command"><pre>${escapeXml(downstream_cmd)}</pre></entry>
</properties>
EOF
"""
}
