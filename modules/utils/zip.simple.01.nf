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


process SIMPLE_ZIP_01 {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${params.prefix?:""}archive.zip"),emit:zip
	path("version.xml"),emit:version
script:
	def prefix = params.prefix?:""
	def compression_level = meta.compression_level?:"5"
"""
hostname 1>&2
set -o pipefail

mkdir -p "${prefix}archive"

cat << EOF | while read F ; do ln -s "\${F}" "./${prefix}archive/" ; done
${L.join("\n")}
EOF

zip -r -${compression_level} "${prefix}archive.zip" "${prefix}archive"

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="id">${task.process}</entry>
	<entry key="num-entries">${L.size()}</entry>
	<entry key="output">${prefix}archive.zip</entry>
</properties>
"""
stub:
"""
touch "${prefix}archive.zip"
echo "<properties/>" > version.xml
"""
}
