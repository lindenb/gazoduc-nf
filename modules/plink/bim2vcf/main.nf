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
process BIM_TO_VCF {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1), path(fasta)
    tuple val(meta2), path(fai)
    tuple val(meta3), path(dict)
    tuple val(meta),path(bim)
output:
    tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.bim2vcf"
    def args  = task.ext.args ?:""
    def awk_expr = task.ext.awk_expr?:"1==1"
    def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
    def jvarkit = task.ext.jvarkit?:"java ${jvm} -jar \${HOME}/jvarkit.jar"
"""
mkdir -p TMP

${jvarkit} dict2vcf "${dict}"  | grep -v '^##FORMAT' > TMP/jeter.vcf

awk -F '\t' '(${awk_expr}) {printf("%s\t%s\t%s\t%s\t%s\t.\t.\t.\\n",\$1,\$4,\$2,\$5,\$6)}' '${bim}' |\\
    ${jvarkit}  bedrenamechr -c 1 -R '${dict}' --convert SKIP >> TMP/jeter.vcf


bcftools norm  --check-ref s --fasta-ref "${fasta}" -O u  TMP/jeter.vcf |\\
    bcftools sort -m ${task.memory.giga}G --temp-dir TMP/sort -O z -o TMP/jeter2.vcf.gz 


bcftools index --threads ${task.cpus} -f -t TMP/jeter2.vcf.gz 

mv TMP/jeter2.vcf.gz  ${prefix}.vcf.gz
mv TMP/jeter2.vcf.gz.tbi  ${prefix}.vcf.gz.tbi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}.bim2vcf"

"""
touch versions.yml "${prefix}.vcf.gz" "${prefix}.vcf.gz.tbi" 
"""
}
