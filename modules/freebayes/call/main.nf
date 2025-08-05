process CALL_FREEBAYES {
tag "${optional_bed?:optiona_bed.name}"
label "single"
afterScript "rm -rf TMP"
array 100
conda "${moduleDir}/../../../conda/freebayes.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path("BAMS/*"),path(optional_bed)
output:
    tuple val(meta),path("*.bcf",arity:'1'),path("*.csi",arity:'1'),path(optional_bed),emit:vcf
    path("versions.yml"),emit:versions
script:
    def args = task.ext.args?:"--use-best-n-alleles 6 "
    def prefix = task.ext.prefix?:"\${MD5}"
"""

mkdir -p TMP
find BAMS/ -name "*am" | sort > TMP/bams.list
test -s  TMP/bams.list

MD5=\$(cat TMP/bams.list ${optional_bed?optional_bed:""} | md5sum | cut -d ' ' -f1)

freebayes \\
        -f "${fasta}" \\
        -L TMP/bams.list \\
        ${optional_bed?"--targets \"${optional_bed}\"":""} \\
        --vcf TMP/jeter.vcf \\
        ${args}
    

awk -F '\t' '/^#/ {print;next;} {OFS="\t";R=\$4;gsub(/[^ATGCatgc]/,"N",R);\$4=R;print;}' TMP/jeter.vcf |\\
bcftools view --threads ${task.cpus} -O b -o TMP/${prefix}.bcf --write-index 
mv TMP/*.bcf ./
mv TMP/*.csi ./


cat << END_VERSIONS > versions.yml
${task.process}:
    freebayes: todo
END_VERSIONS
"""
}
