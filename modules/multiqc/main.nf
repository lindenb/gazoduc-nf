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
process MULTIQC {
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/multiqc.yml"
input:
	tuple val(meta),path(multiqc_files, stageAs: "?/*")
output:
	tuple val(meta),path("*.zip"),optional:true,emit:zip
	tuple val(meta),path("*multiqc/*_data"),optional:true,emit:datadir
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:""
    def comment = task.ext.comment ?:""
    def title = task.ext.comment ?:""
"""
mkdir -p TMP
export TMPDIR=\${PWD}/TMP

mkdir -p "${prefix}multiqc"

export LC_ALL=en_US.utf8

multiqc --filename  "${prefix}multiqc_report.html" \\
    --no-ansi \\
	--title "${title}" \\
	--comment "${comment}"  \\
	--force \\
	--outdir "${prefix}multiqc"  \\
	.
		
zip -9 -r "${prefix}multiqc.zip" "${prefix}multiqc"


cat << EOF > versions.yml
${task.process}:
    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
EOF
"""

stub:
"""
mkdir -p "${meta.id}multiqc/${meta.id}_data"
touch "${meta.id}multiqc/${meta.id}_data/jeter.json"
touch multiqc.zip
touch versions.yml
"""
}