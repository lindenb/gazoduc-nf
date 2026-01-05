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
process TABIX_BED {
    label "process_single"
	tag "${meta.id}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    input:
		tuple val(meta ),path(bed)
    output:
		tuple val(meta),path("*.gz"),path("*.gz.tbi"),emit: bed
		path("versions.yml"),emit:versions
    script:
		def prefix= task.ext.prefix?:(bed.name.endsWith(".bed")?bed:bed.baseName)
    """
	hostname 1>&2
	mkdir -p TMP
	${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" |\\
        C_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k1,1 -k2,2n -T TMP |\\
        bgzip > ${prefix}.gz

	tabix -p bed -f ${prefix}.gz


cat << END_VERSIONS > versions.yml
"${task.process}":
	tabix: \$(tabix  --version |  awk '(NR==1){print \$NF}')
END_VERSIONS
    """

stub:
def prefix= task.ext.prefix?:(bed.name.endsWith(".bed")?bed:bed.baseName)
"""
touch versions.yml ${prefix}.gz ${prefix}.gz.tbi
"""
}
