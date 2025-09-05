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
/** ADD arbitrary Genotype to CandidateSV VCF */
process MANTA_CANDIDATESV {
    label "process_single"
    tag "${meta.id?:""} ${vcf.name}"
    afterScript "rm -rf TMP"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    input:
        tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    	path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:"${meta.id}.GT.candidateSV"
	"""
	hostname 1>&2
    bcftools view "${vcf}" |\\
        awk -F '\t' '/^##/ {print;next;} /^#CHROM/ {printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"original manta ${vcf} doesn not contain a genotype. This is an arbitratry Genotype for downstream analysis. All Genotypes will be HET\\">\\n%s\tFORMAT\t${meta.id}\\n",\$0);next;} {printf("%s\tGT\t0/1\\n",\$0);}' |\\
        bcftools view -O z -o "${prefix}.vcf.gz"
    
    bcftools index -f -t "${prefix}.vcf.gz"

cat << EOF > versions.yml
"${task.process}":
	bcftools: TODO
EOF
"""
}


