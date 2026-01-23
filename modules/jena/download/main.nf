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

process JENA_DOWNLOAD {
label "process_single"
input:
	val(meta)
output:
	tuple val(meta),path("apache-jena/bin/arq"),emit:arq
    tuple val(meta),path("apache-jena/bin/riot"),emit:riot
    tuple val(meta),path("apache-jena/bin/tdbquery"),emit:tdbquery
    tuple val(meta),path("apache-jena/bin/tdbloader"),emit:tdbloader
	path("versions.yml"),emit:versions
script:
    def version = task.ext.version?:"5.6.0"
	def url = "https://dlcdn.apache.org/jena/binaries/apache-jena-${version}.tar.gz"
"""
curl -L -o apache-jena-${version}.tar.gz "${url}"
tar xvfz apache-jena-${version}.tar.gz

mv -v apache-jena-${version} apache-jena
find apache-jena 1>&2

cat << EOF > versions.yml
"${task.process}"
	jena: "${url}"
EOF
"""

stub:
"""
touch versions.yml
mkdir -p apache-jena/bin
(cd apache-jena/bin && touch arq riot tdbquery tdbloader)
"""
}
