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
process CIRCOS {
tag "${fasta.name}"
label "process_medium"
conda "${moduleDir}/../../../conda/circos.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta ),path(config)
    tuple val(meta2),path("DATA/*")
output:
	tuple val(meta),path("*.svg"),optional:true,emit:svg
    tuple val(meta),path("*.png"),optional:true,emit:png
script:
    def prefix = task.ext.prefix ?: "${meta.id}."
    def args1 = task.ext.args1?:""
"""
mkdir -p TMP
circos --conf ${config} ${args1} -outputdir TMP

mv TMP/circos.svg "${prefix}circos.svg"
mv TMP/circos.png "${prefix}circos.png"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    circos: "\$(circos -v 2>&1)"
END_VERSIONS
"""

stub:
"""
touch versions.yml "${meta.id}.png" "${meta.id}.svg"
"""
}
