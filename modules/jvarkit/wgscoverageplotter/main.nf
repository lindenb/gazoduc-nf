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

process WGS_COVERAGE_PLOTTER  {
    tag "${meta.id}"
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	input:
		tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
		tuple val(meta ),path(bam),path(bai)
	output:
		tuple val(meta ),path("*.{svg,svg.gz}"),emit:svg
		path("versions.yml"),emit:versions
	script:
		def args = task.ext.args?:""
        def prefix = task.ext.prefix?:"${meta.id}"
		def jvm = " -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
set -o pipefail

mkdir -p TMP
jvarkit ${jvm} wgscoverageplotter \\
    ${args} \\
    -R "${fasta}" "${bam}" > TMP/jeter.svg
    
mv -v TMP/jeter.svg "${prefix}.svg"

cat << EOF > versions.yml
${task.process}:
	jvarkit: "\$(jvarkit --version)"
EOF
"""

stub:
   def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml  ${prefix}.svg
"""

}
