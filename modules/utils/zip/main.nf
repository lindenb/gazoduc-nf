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


process ZIP {
label "process_single"
tag "${meta.id}"
input:
	tuple val(meta),path(L)
output:
	tuple val(meta),path("*.zip"),emit:zip
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
	def compression_level = task.ext.level?:"5"
"""


mkdir -p "${prefix}.archive"

cat << EOF | while read F ; do ln -v -s "\${PWD}/\${F}" "./${prefix}.archive/" ; done
${L.join("\n")}
EOF

zip -r -${compression_level} "${prefix}.archive.zip" "${prefix}.archive/"

cat << EOF > versions.yml
"${task.process}"
	zip: \$(zip --version | grep "This is")
EOF
"""
stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch "${prefix}.archive.zip" versions.yml
"""
}
