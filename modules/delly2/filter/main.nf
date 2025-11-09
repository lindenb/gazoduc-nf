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
process FILTER_DELLY {
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/delly.yml"
    label "process_short"
    input:
		tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
		tuple val(meta),path(vcfin),path(bci)
    output:
		tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
        path("*.md5"),emit:md5
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:"${meta.id}.sv"
    """
    export LC_ALL=C
    export PATH=\${PWD}:\${PATH}
    mkdir -p TMP

    delly filter ${tag} -f germline -o TMP/jeter.bcf "${vcfin}" 1>&2

    bcftools sort --max-mem "${task.memory.giga}G" -T TMP/sort -O v -o "TMP/jeter1.vcf" TMP/jeter.bcf

    bcftools  +fill-tags -O v  -o TMP/jeter2.vcf TMP/jeter1.vcf  -- -t AN,AC,AF,AC_Hom,AC_Het,AC_Hemi,NS
    mv TMP/jeter2.vcf TMP/jeter1.vcf

    bcftools view -O b9 -o "TMP/${prefix}.bcf" TMP/jeter1.vcf
    bcftools index -f "TMP/${prefix}.bcf"

    mv TMP/${prefix}.* ./

    md5sum ${prefix}.bcf >  ${prefix}.bcf.md5

    touch versions.yml
    """

    stub:
     def prefix=task.ext.prefix?:"${meta.id}.sv"
    """
    touch versions.yml ${prefix}.bcf ${prefix}.bcf.csi ${prefix}.bcf.md5
    """
    }
