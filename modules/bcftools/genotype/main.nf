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
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



*/
process BCFTOOLS_GENOTYPE {
    label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val( meta1),path(fasta)
        tuple val(_meta2),path(_fai)
        tuple val( meta3),path(alleles),path(alleles_tbi) //CHROM(tab)POS(tab)REF,ALT cf. bcftools manual
        tuple val( meta ),path("BAMS/*")
output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
        path("versions.yml"),emit:versions
script:
        def prefix = task.ext.prefix?:"${meta.id}.${meta3.id}"
        def mapq = task.ext.mapq?:10
        def ploidy = "";
        if(meta1.ucsc_name=="hg19") {
                ploidy="GRCh37"
                }
        else if(meta1.ucsc_name=="hg38") {
                ploidy="GRCh38"
                }
        def args1 = task.ext.args1?:""
        def args2 = task.ext.args2?:""
        def args3 = task.ext.args3?:""
	def annotations1 = task.ext.annotations1?:"-a 'FORMAT/AD' -a 'FORMAT/DP' -a 'FORMAT/SP' -a 'FORMAT/SCR' -a 'FORMAT/QS' -e 'FORMAT/NMBZ'"
	def annotations2 = task.ext.annotations2?:"-a 'INFO/PV4' -a 'FORMAT/GQ' -a 'FORMAT/GP'"
"""
hostname 1>&2
mkdir -p TMP
find BAMS -name "*.bam" -o -name "*.cram" |\\
        samtools samples |\\
        sort -t '\t' -T TMP -k1,1 |\\
        cut -f2 > TMP/jeter.list

test -s  TMP/jeter.list

bcftools mpileup \\
    ${args2} \\
    --redo-BAQ \\
    ${annotations1} \\
    --threads ${task.cpus} \\
    --fasta-ref "${fasta}" \\
    --regions-file "${alleles}" \\
    -q ${mapq} \\
    -O u \\
    -o TMP/jeter2.bcf \\
    --bam-list TMP/jeter.list

bcftools call  \\
    ${args3} \\
    --keep-alts \\
    ${annotations2} \\
    ${ploidy.isEmpty()?"":"--ploidy ${ploidy}"} \\
    --threads ${task.cpus}  \\
    --targets-file "${alleles}" \\
    --constrain alleles \\
    --multiallelic-caller \\
    --output-type z -o "TMP/jeter3.vcf.gz" TMP/jeter2.bcf

mv -v TMP/jeter3.vcf.gz "${prefix}.vcf.gz"

bcftools index --threads ${task.cpus} -t -f "${prefix}.vcf.gz"

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""
stub:

def prefix = task.ext.prefix?:"${meta.id}.${meta3.id}"
"""
touch versions.yml "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi"
"""
}
