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
 * genotype one or more bam for variants in a VCF
 */
process GATK_GENOTYPE_VCF {
label "process_single"
tag "bams:${meta.id} vcf:${metaV.id}"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
    tuple val(meta4),path(dbsnp),path(dbsnp_tbi) //optional
    tuple val(metaV),path(vcf),path(vcf_id)
	tuple val(meta ),path("BAMS/*")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
	def prefix = task.ext.meta?:"${meta.id}.${metaV.id}"
    def args1 = task.ext.args1?:""
    def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true "
"""
mkdir -p TMP
find BAMS/ -name "*.bam" -o -name "*.cram" > TMP/jeter.list
test -s TMP/jeter.list

gatk --java-options "${jvm}" HaplotypeCaller \\
    ${args1} \\
	-I TMP/jeter.list \\
	-L ${vcf} \\
    ${dbsnp?"--dbsnp \"${dbsnp}\"":""} \\
	-R "${fasta}" \
	--alleles ${vcf}\\
	--output-mode EMIT_ALL_CONFIDENT_SITES \\
	-O "TMP/jeter3.vcf.gz"

# remove 'Infinity' from QUAL
gunzip -c TMP/jeter3.vcf.gz |\\
	awk -F '\t' '/^#/ {print; OFS="\t";next;} {if(\$6=="Infinity") {\$6=".";} print;}' |\\
	bcftools view -O z -o TMP/jeter4.vcf.gz

bcftools index -ft TMP/jeter4.vcf.gz

mv TMP/jeter4.vcf.gz ./${prefix}.vcf.gz
mv TMP/jeter4.vcf.gz.tbi  ./${prefix}.vcf.gz.tbi

cat << EOF > versions.yml
${task.process}:
    gatk: "\$((gatk --java-options "${jvm}" --version 2> /dev/null  | paste -s -d ' ' ) || true)"
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
EOF
"""
stub:
	def prefix = task.ext.meta?:"${meta.id}.${metaV.id}"
"""
touch versions.yml "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi"
"""
}
