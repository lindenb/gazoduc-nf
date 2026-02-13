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
process PLINK_BFILE2VCF {
tag "${meta.id}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
//afterScript "rm -rf TMP"
input:
    tuple val(meta1), path(fasta) //optional
    tuple val(meta2), path(fai) // optional
    tuple val(meta2), path(dict) // optional
    tuple val(meta),path(bim),path(bed),path(fam)
output:
    tuple val(meta),path("*.vcf.gz"),path("*.tbi"),emit:vcf
    tuple val(meta),path("*.log"),emit:log
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def args1 = task.ext.args1?:" --allow-extra-chr --allow-no-sex  --threads ${task.cpus}  --const-fid 1  --memory ${task.memory.mega}"
    def args2 = task.ext.args2?:"--keep-allele-order"
    def jvm = task.ext.jvm?:"-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
mkdir TMP

plink2 \\
    ${fasta ?"--ref-from-fa force --fa ${fasta} --real-ref-alleles ":""} \\
    ${args1} \\
    --bim ${bim} \\
    --bed ${bed} \\
    --fam ${fam} \\
    --export vcf bgz id-delim='#' \\
    ${args2} \\
    --out TMP/jeter

#
# rename samples, remove '1#' in front of VCF
#
bcftools query -l TMP/jeter.vcf.gz > TMP/jeter.a
sed 's/^1#//' TMP/jeter.a > TMP/jeter.b
sort -T TMP TMP/jeter.b | uniq -d > TMP/jeter.dups
test ! -s TMP/jeter.dups
paste TMP/jeter.a TMP/jeter.b > TMP/rename.tsv

bcftools reheader  --threads ${task.cpus}  --temp-prefix TMP/rh --samples TMP/rename.tsv -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
 mv -v TMP/jeter2.vcf.gz TMP/jeter.vcf.gz

if ${fasta?true:false}
then


   
    gunzip -c TMP/jeter.vcf.gz |\\
    awk -F '\t' '
               BEGIN{OFS="\t";}
               /^#/ {print;next;}
               (\$4=="0") {next;}
               (\$5=="0") {next;}
               (\$4=="-") {\$2=sprintf("%d",int(\$2)-1);\$4="N";\$5=sprintf("N%s",\$5); print;next;}
               (\$5=="-") {\$2=sprintf("%d",int(\$2)-1);\$5="N";\$4=sprintf("N%s",\$4); print;next;}
                          {print;}' |\\
    jvarkit  ${jvm} vcfsetdict \\
        -n SKIP \\
       -R "${fasta}" |\\
    bcftools norm --fasta-ref ${fasta} --check-ref s -O z -o TMP/jeter2.vcf.gz

    mv -v TMP/jeter2.vcf.gz TMP/jeter.vcf.gz


fi

bcftools sort --max-mem '${task.memory.giga}G' -T TMP/tmp -O b -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz

bcftools index --threads ${task.cpus}  --force --tbi TMP/jeter2.vcf.gz


mv -v TMP/jeter.log "${meta.id}.log"
mv -v TMP/jeter2.vcf.gz ${meta.id}.vcf.gz
mv -v TMP/jeter2.vcf.gz.tbi ${meta.id}.vcf.gz.tbi

touch versions.yml
"""

stub:
"""
touch versions.yml "${meta.id}.vcf.gz" "${meta.id}.vcf.gz.tbi" ${meta.id}.log"
"""
}
