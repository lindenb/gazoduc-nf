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
include {verify } from '../../../modules/utils/functions.nf'

process VCF_STRECH_TO_SVG  {
    tag "${meta.id}"
	label "process_single"
	conda "${moduleDir}/../../../conda/bioinfo.01.yml"
	afterScript "rm -rf TMP"
	input:
		tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
		tuple val(meta4),path(gtf),path(gtf_tbi)
		tuple val(meta5 ),path("BAMS/*")
		tuple val(meta6  ),path(vcf),path(tbi)
		tuple val(meta),path(bed)
	output:
		tuple val(meta ),path("*.{svg,svg.gz}"),emit:svg
		path("versions.yml"),emit:versions
	script:
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:"--hom-ref"
        //def prefix = task.ext.prefix?:"${meta.id}"
		def jvm = "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
set -o pipefail

mkdir -p TMP/OUT BAMS
find BAMS/ -\\( -name "*.bam" -o -name "*.cram" \\) | sort -V -T TMP > TMP/bams.list

jvarkit ${jvm} vcfstrech2svg \\
    ${args1} \\
	${args2} \\
	--bed "${bed}" \\
	--bam-list TMP/bams.list \\
	${gtf?"--gtf \"${gtf}\" ":""} \\
    -R "${fasta}" \\
	-o TMP/OUT \\
	"${vcf}"
    
mv -v TMP/OUT/*.svg ./ || true
mv -v TMP/OUT/*.svg.gz ./ || true

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
