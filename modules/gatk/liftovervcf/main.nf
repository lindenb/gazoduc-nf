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

process LIFTOVER_VCF {
tag "${meta.id?:vcf.name} ${chain.name}"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
label "process_short"
afterScript 'rm -rf TMP'
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta4),path(chain)
	tuple val(meta ),path(vcf),path(idx)
output:
	tuple val(meta),path("*.lift.vcf.gz"),path("*.lift.vcf.gz.tbi"),emit:vcf
    tuple val(meta),path("*.fail.vcf.gz"),path("*.fail.vcf.gz.tbi"),emit:fail
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
    def prefix2 = "${prefix}${meta1.ucsc_name?".${meta1.ucsc_name}":""}"
	def args1 =  task.ext.args1?:""
    def args2 =  task.ext.args1?:""
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP" 
"""
hostname 1>&2
mkdir -p TMP

if ${vcf.name.endsWith(".bcf")}
then
    bcftools view  ${args1} --threads ${task.cpus} -Oz -o TMP/jeter1.vcf.gz "${vcf}" 
fi

gatk --java-options "${jvm}" LiftoverVcf \\
	--CHAIN ${chain} \\
	--INPUT ${vcf.name.endsWith(".bcf")?"TMP/jeter1.vcf.gz":"${vcf}"} \\
	--OUTPUT TMP/jeter2.vcf.gz \\
	--REFERENCE_SEQUENCE ${fasta} \\
	--REJECT TMP/reject.vcf.gz \\
    --DISABLE_SORT true \\
    ${args2}


bcftools sort --max-mem '${task.memory.giga}G' \\
    -T TMP/tmp -O z -o TMP/jeter3.vcf.gz TMP/jeter2.vcf.gz
bcftools index --threads '${task.cpus}' -t --force TMP/jeter3.vcf.gz

bcftools sort --max-mem '${task.memory.giga}G' \\
    -T TMP/tmp -O z -o TMP/jeter4.vcf.gz TMP/reject.vcf.gz
bcftools index --threads '${task.cpus}' -t --force TMP/jeter4.vcf.gz

mv TMP/jeter3.vcf.gz "${prefix2}.lift.vcf.gz"
mv TMP/jeter3.vcf.gz.tbi "${prefix2}.lift.vcf.gz.tbi"

mv TMP/jeter4.vcf.gz "${prefix}.fail.vcf.gz"
mv TMP/jeter4.vcf.gz.tbi "${prefix}.fail.vcf.gz.tbi"
 
cat << EOF > versions.yml
${task.process}:
    gatk: "\$( gatk --version 2> /dev/null  | paste -s -d ' ' )"
EOF
"""

stub:
"""
touch "${meta.id}.lift.vcf.gz"
touch "${meta.id}.lift.vcf.gz.tbi"
touch "${meta.id}.fail.vcf.gz"
touch "${meta.id}.fail.vcf.gz.tbi"
touch versions.yml
"""
}
