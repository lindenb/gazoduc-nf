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


process SORT_VCF {
label "process_short"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
    tuple val(meta ),path(vcf)
output:
    tuple val(meta),path("*.vcf.gz"), path("*.vcf.gz.tbi"),emit:vcf
    tuple val(meta),path("*.md5"),optional:true,emit:md5
    path("versions.yml"),emit:versions
script:
   def prefix = task.ext.prefix?:"${meta.id}.sort"
   def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
   def with_md5 = task.ext.with_md5?:true
   def compression_level = task.ext.compression_level?:"9"
"""	
	hostname 1>&2
	mkdir -p TMP

    gatk --java-options "${jvm}" SortVcf \\
        --INPUT ${vcf} \\
        --COMPRESSION_LEVEL ${compression_level} \\
        --OUTPUT TMP/jeter.vcf.gz \\
        --TMP_DIR \${PWD}/TMP \\
        --CREATE_INDEX false
    
    bcftools index --threads ${task.cpus} -f -t TMP/jeter.vcf.gz

    mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
    mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi 

    # Generate MD5 if needed
    if ${with_md5}
    then 
        md5sum ${prefix}.vcf.gz > ${prefix}.vcf.gz.md5
    fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""

stub:
   def prefix = task.ext.prefix?:"${meta.id}.sort"
"""
touch versions.yml ${prefix}.vcf.gz  ${prefix}.vcf.gz.tbi   ${prefix}.vcf.gz.md5
"""
}
