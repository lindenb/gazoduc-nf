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
process MERGE_GENOTYPES {
    tag "${meta.id}"
    label "process_short"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/delly.yml"
    input:
	    tuple val(meta),path("VCFS/*")
    output:
	    tuple val(meta),path("*.bcf",arity:"1"),path("*.csi",arity:"1"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:(meta.id?:"merged")
    // note to self , cd gazoduc-nf if not enough memory
    """
    hostname 1>&2
    ulimit -s unlimited || true
    mkdir -p TMP
    find VCFS/ -name "*.bcf" | sort > TMP/jeter.list
    test -s TMP/jeter.list

    bcftools merge --force-single --threads ${task.cpus} -m id -O b9 -o TMP/${prefix}.bcf --file-list TMP/jeter.list
    bcftools index --threads ${task.cpus} --csi TMP/${prefix}.bcf 
    
    mv TMP/${prefix}.bcf ./
    mv TMP/${prefix}.bcf.csi ./

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1) {print \$NF;}')
END_VERSIONS
    """

stub:
 def prefix = task.ext.prefix?:(meta.id?:"merged")
"""
touch versions.yml {prefix}.bcf {prefix}.bcf.csi
"""
}