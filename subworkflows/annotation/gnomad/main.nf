
/**


UNDER CONSTRUCTION



*/
workflow GNOMAD {
meta:
    meta
    fasta
    fai
    dict
    vcf
main:
    ANNOT(vcf)
emit:
    vcf = ANNOT.out.vcf
}

process ANNOT {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"

input:
    tuple val(meta),path(vcf),path(idx)
output:
    tuple val(meta),path("*.bcf"),path("*.csi"),emit:vcf
script:
"""
bcftools view "${vcf}" |\\
jvarkit -Xmx${task.memory.giga}G -XX:-UsePerfData  -Djava.io.tmpdir=TMP vcfgnomad  \\
    --bufferSize 10000 \
    --gnomad "${gnomadVcf}" \
    --fields "${gnomadPop}" \
    --max-af "${gnomadAF}" > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf

## bcftools ne marche pas avec   -e 'FILTER ~ "GNOMAD_GENOME_BAD_AF" || FILTER ~ "GNOMAD_GENOME_AS_VQSR"' filtre pas forcement 
awk  '\$0 ~ /^#/ || !(\$7 ~ /GNOMAD_GENOME_BAD/)' TMP/jeter1.vcf > TMP/jeter2.vcf
mv TMP/jeter2.vcf TMP/jeter1.vcf
"""
}