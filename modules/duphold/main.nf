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
process DUPHOLD {
tag "${meta.id?:""}"
label 'process_single'
conda "${moduleDir}/../../conda/duphold.yml"
afterScript "rm -rf TMP"

input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta), path(bam), path(bai), path(svvcf), path(optional_snp) /* bug in duphold with optional_vcf */, path(optional_snp_idx)
output:
    tuple val(meta), path("*.vcf.gz") ,path("*.tbi"),emit: vcf
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
"""
hostname 1>&2
mkdir -p TMP


bcftools view --header-only ${svvcf} > TMP/vcf.header
bcftools query -l  TMP/vcf.header | sort -T TMP | uniq > TMP/vcf.samples
comm -12 TMP/vcf.samples <(echo "${meta.id}") > TMP/vcf.contains.sample.flag

if test -s vcf.contains.sample.flag
then
    bcftools view --threads ${task.cpus} -O b --write-index -o TMP/jeter.bcf "${svvcf}"
else
    bcftools view -G -O v "${svvcf}" |\\
        awk -F '\t' '/^##/ {print;next;} /^#CHROM/ {printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\\n");printf("%s\tFORMAT\t${meta.id}\\n",\$0);next;} {printf("%s\tGT\t0/1\\n",\$0);}' |\\
        bcftools view -O b -o --write-index -o TMP/jeter.bcf
fi


duphold \\
    ${args} \\
    ${optional_snp ? "--snp \"${optional_snp}\"" : ""} \\
    --threads ${task.cpus} \\
    --output TMP/jeter.vcf.gz \\
    --vcf TMP/jeter.bcf \\
    --bam "${bam}" \\
    --fasta "${fasta}"

bcftools index -f -t --threads "${task.cpus}" TMP/jeter.vcf.gz

mv  TMP/jeter.vcf.gz ${prefix}.vcf.gz
mv  TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    duphold: \$(duphold -h | head -n 1 | sed -e "s/^version: //")
END_VERSIONS
"""
}
