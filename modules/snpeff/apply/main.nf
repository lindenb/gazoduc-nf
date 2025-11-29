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
process SNPEFF_APPLY {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"

input:
	tuple val(meta1 ),path(snpeff_datadir),path(snpeff_db_name)
	tuple val(meta ),path(vcf)
output:
	tuple val(meta),path("*.vcf.gz"),optional:true,emit:vcf
	path("versions.yml"),emit:versions
script:
	def args1=task.ext.args1?:""
	def args2=task.ext.args2?:"-nodownload"
	def args3=task.ext.args3?:""
	def jvm = task.ext.jvm?:" -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
	def prefix = task.ext.prefix?:"${meta.id}"
	def accessions  = task.ext.accessions?:""
"""
mkdir -p TMP

bcftools view ${args1} "${vcf}" |\\
	snpEff ${jvm}  eff \\
        -dataDir "\${PWD}/${snpeff_datadir}" \\
		${args2} \\
        `cat ${snpeff_db_name}` |\\
	bcftools view -O z -o TMP/jeter.vcf.gz


if ${!accessions.trim().isEmpty()}
then
   bcftools view TMP/jeter.vcf.gz |\\
   jvarkit ${jvm}  vcffilterso \\
        ${args3} \\
        --acn "${accessions}" |\\
	bcftools view -O z -o TMP/jeter2.vcf.gz
    
    mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
fi

if test \$(bcftools query -f "\\n" TMP/jeter.vcf.gz |wc -l) -gt 0
then
	mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
fi

##################
cat << EOF > versions.yml
"${task.process}":
	snpeff: todo
	jvarkit: todo
EOF
"""

stub:
	def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.vcf.gz
"""
}
