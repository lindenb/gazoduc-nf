/*

Copyright (c) 2024 Pierre Lindenbaum

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
process HC_GENOTYPE {
tag "${bed.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
	tuple val(meta4),path(dbsnp)
	tuple val(meta5),path(dbsnp_tbi)
        tuple path(bed),path(gvcf),path(gvcf_idx)
output:
	tuple path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:output
script:
"""
hostname 1>&2
mkdir -p TMP



	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \\
		GenotypeGVCFs \\
	        ${dbsnp?"--dbsnp ${dbsnp}":""} \\
		-R "${fasta}" \\
		-L "${bed}" \\
		-O TMP/genotyped.vcf.gz \\
		-V "${gvcf}" \\
		-G StandardAnnotation \\
		-G AS_StandardAnnotation



mv -v TMP/genotyped.vcf.gz "${bed.getBaseName()}.vcf.gz"
mv -v TMP/genotyped.vcf.gz.tbi "${bed.getBaseName()}.vcf.gz.tbi"
"""
}

