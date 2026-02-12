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
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



*/

/**
 * emit vcf with message EMPTY or VARIANT
 */
process BCFTOOLS_COUNTS {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta),path(vcf)
    output:
        tuple val(meta),path("*.tsv"),emit:counts
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:"${meta.id}"
        def args = task.ext.args?:""
    """
    mkdir -p TMP
    bcftools +counts ${args}  ${vcf} |\\
        sed 's/^Number of //;s/:[ ]*/\t/' > TMP/jeter.tsv

    echo -n "id\tvcf\t" > TMP/jeter2.tsv
    cut -f 1 TMP/jeter.tsv | paste -sd '\t'  >> TMP/jeter2.tsv

    echo -n "${meta.id}\t${vcf.name}\t" >> TMP/jeter2.tsv
    cut -f 2 TMP/jeter.tsv | paste -sd '\t'  >> TMP/jeter2.tsv

    mv TMP/jeter2.tsv "${prefix}.tsv"

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """

stub:
     def prefix = task.ext.prefix?:"${meta.id}"
    """
    touch versions.yml "${prefix}.tsv"
    """
}
