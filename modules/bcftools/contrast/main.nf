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
/**
 * apply bctools contrast 
 */
process BCFTOOLS_CONTRAST {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1),path(cases)
        tuple val(meta2),path(controls)
        tuple val(meta),path(vcf)
    output:
        tuple val(meta),path("*.vcf.gz"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:"${meta.id}.contrast"
        def annotations = task.ext.annotations?:"PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT"
        def args1 = task.ext.args1?:""
    """
    mkdir -p TMP
    
	bcftools +contrast \\
		--case-samples ./${cases} \\
		--control-samples ./${controls} \\
        ${args1} \\
        -a "${annotations}" \\
		-O z -o TMP/jeter.vcf.gz '${vcf}'

mv  TMP/jeter.vcf.gz ${prefix}.vcf.gz


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
stub:
        def prefix = task.ext.prefix?:"${meta.id}.contrast"
"""
touch versions.yml ${prefix}.vcf.gz
"""
}
