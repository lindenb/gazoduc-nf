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
process GHOSTSCRIPT_MERGE {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/ghostscript.yml"
afterScript "rm -rf TMP"
input:
        tuple val(meta),path("PDF/*")
output:
        tuple val(meta),path("*.pdf"),emit:pdf
        path("versions.yml"),emit:versions
script:
        def cmd = task.ext.cmd?:"| sort -V -T TMP"
        def prefix  = task.ext.prefix?:meta.id
"""
hostname 1>&2
mkdir -p TMP

find PDF/ -type l -name "*.pdf" ${cmd} > TMP/jeter.txt

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=TMP/jeter.pdf @TMP/jeter.txt

mv -v TMP/jeter.pdf "${prefix}.pdf"

cat << END_VERSIONS > versions.yml
"${task.process}":
	gs: "\$(gs --version )"
END_VERSIONS
"""

stub:
"""
touch "${meta.id}.pdf" versions.yml
"""
}
