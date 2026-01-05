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
process XHUNTER_DENOVO_MERGE {
    label "process_single"
	tag "${meta3.id}"
	afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/xhunter.denovo.yml"
	input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
		tuple val(meta3),path("JSON/*")
        tuple val(meta4),path(manifest)
	output:
		tuple val(meta3),path("*.json"),emit:json
		//tuple val(meta3),path("manifest.txt"),emit:manifest
		path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta3.id}"
"""
hostname 1>&2
mkdir -p TMP

ExpansionHunterDenovo merge \\
		--manifest ${manifest} \\
		--reference "${fasta}" \\
		--output-prefix "TMP/${prefix}" 1>&2

mv -v TMP/*.json ./

cat << EOF > versions.yml
"${task.process}":
    xhunter: \$(ExpansionHunterDenovo --version 2>&1 | awk '{print \$3}')
EOF
"""

stub:
"""
touch versions.yml ${meta3.id}.json
"""
}
