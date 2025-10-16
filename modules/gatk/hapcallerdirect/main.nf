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

process HAPLOTYPECALLER_DIRECT {
tag "${meta.id?:""} ${optional_bed?optional_bed.name:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(dbsnp),path(dbsnptbi)
    tuple val(meta5),path(pedigree)
    tuple val(meta ),path("BAMS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),path(optional_bed),emit:vcf
    path("versions.yml"),emit:versions
script:
   def prefix0 = (meta.id?:"\${MD5}")
   def prefix= task.ext.prefix?:prefix0
   def args1 = task.ext.args1?:"-G StandardAnnotation -G StandardHCAnnotation"
   def args2 = task.ext.args2?:""
   def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
"""
hostname 1>&2
mkdir -p TMP
set -x
find ./BAMS/ \\( -name "*.bam" -o -name "*.cram" \\) | sort -T TMP -V >> TMP/bams.list
MD5=\$(cat TMP/bams.list ${optional_bed?optional_bed:""} | md5sum | cut -d ' ' -f1)

cat TMP/bams.list | samtools samples -f "${fasta}" | awk -F '\t' '(\$3==".")' > TMP/bad.txt
if test -s TMP/bad.txt
then
    sed 's/^/UNDEFINED or WRONG REF FOR /' TMP/bad.txt 1>&2
fi

test ! -s TMP/bad.txt

gatk --java-options "${jvm}" HaplotypeCaller \\
    ${optional_bed?"-L ${optional_bed}":""} \\
    ${dbsnp?"--dbsnp \"${dbsnp}\"":""} \\
    ${pedigree?"--pedigree \"${pedigree}\"":""} \\
    -R "${fasta}" \\
    -I TMP/bams.list \\
    ${args1} ${args2}\\
    -O "TMP/jeter.vcf.gz"


bcftools index --threads ${task.cpus} --force --tbi "TMP/jeter.vcf.gz"

mv "TMP/jeter.vcf.gz" ${prefix}.vcf.gz
mv "TMP/jeter.vcf.gz.tbi" ${prefix}.vcf.gz.tbi

cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""

stub:
 def prefix = (meta.id?:(optional_bed?optional_bed.name:""))
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
}
