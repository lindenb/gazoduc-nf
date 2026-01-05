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


process FOP_DOWNLOAD {
tag "2.11"
executor "single"
input:
	val(meta)
output:
	tuple val(meta),path("fop-2.11/fop/fop"),emit:fop
	path("versions.yml"),emit:versions
script:
	def url="https://www.apache.org/dyn/closer.cgi?filename=/xmlgraphics/fop/binaries/fop-2.11-bin.tar.gz&action=download
"""
curl -L -o jeter.zip "${url}"
unzip jeter.zip
rm jeter.zip

cat << EOF > versions.yml
"${task.process}"
	fop: "${url}"
EOF
"""
}
