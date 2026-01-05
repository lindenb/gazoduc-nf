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
process CALL_FREEBAYES {
tag "${optional_bed?:optional_bed.name}"
label "process_single"
label "array100"
afterScript "rm -rf TMP"

conda "${moduleDir}/../../../conda/freebayes.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path("BAMS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.bcf",arity:'1'),path("*.csi",arity:'1'),path(optional_bed),emit:vcf
    path("versions.yml"),emit:versions
script:
    def args = task.ext.args?:"--use-best-n-alleles 6 --limit-coverage 500"
    def prefix = task.ext.prefix?:"\${MD5}"
"""

mkdir -p TMP
find BAMS/ -name "*am" | sort > TMP/bams.list
test -s  TMP/bams.list

MD5=\$(cat TMP/bams.list ${optional_bed?optional_bed:""} | md5sum | cut -d ' ' -f1)

freebayes \\
        -f "${fasta}" \\
        -L TMP/bams.list \\
        --genotype-qualities \\
        ${optional_bed?"--targets \"${optional_bed}\"":""} \\
        --vcf TMP/jeter.vcf \\
        ${args}
    

awk -F '\t' '/^#/ {print;next;} {OFS="\t";R=\$4;gsub(/[RYMKSWHBVD]/,"N",R);\$4=R; R=\$5;gsub(/[RYMKSWHBVD]/,"N",R);\$5=R;print;}' TMP/jeter.vcf |\\
bcftools sort --max-mem ${task.memory.giga}G  -T TMP/sort -O b -o TMP/${prefix}.bcf

bcftools index -f --threads ${task.cpus} TMP/${prefix}.bcf

mv TMP/*.bcf ./
mv TMP/*.csi ./


cat << END_VERSIONS > versions.yml
${task.process}:
    freebayes: todo
END_VERSIONS
"""

stub:
"""
MD5=\$(find BAMS/ && echo  ${optional_bed?optional_bed:""} | md5sum | cut -d ' ' -f1)
touch versions.yml "\${MD5}.bcf" "\${MD5}.bcf.csi" 
"""
}
