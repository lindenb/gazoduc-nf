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

process COMBINEGVCFS {
tag "${meta.id?:""} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta),path("VCFS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.g.vcf.gz"),path("*.g.vcf.gz.tbi"),path(optional_bed),emit:gvcf
    path("versions.yml"),emit:versions
script:
   def prefix= task.ext.prefix?:"\${MD5}"
   def args1 = task.ext.args1?:"-G StandardAnnotation -G AS_StandardAnnotation"
"""
hostname 1>&2
mkdir -p TMP

find VCFS/ -name "*.g.vcf.gz" | sort > TMP/jeter.list
MD5=`cat TMP/jeter.list ${optional_bed?:""} | md5sum | awk '{print \$1}'`

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
    CombineGVCFs \\
    -R "${fasta}" \\
    ${optional_bed?"-L \"${optional_bed}\"":""} \\
    -V TMP/jeter.list \\
    -O "TMP/combine.g.vcf.gz" \\
    ${args1}

mv TMP/combine.g.vcf.gz "${prefix}.g.vcf.gz"
mv TMP/combine.g.vcf.gz.tbi "${prefix}.g.vcf.gz.tbi"

cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""
}
