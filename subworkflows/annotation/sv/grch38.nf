
workflow VEP {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        versions = Channel.empty()
        APPLY_VEP(
            fasta,
            fai,
            vcf
            )
    emit:
        vcf=APPLY_VEP.out.vcf
        versions
}

process APPLY_VEP {
tag "${meta.id} ${vcf.baseName}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path(vcf),path(tbi)
output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
script:
    def assembly = task.ext.assembly?:"GRCh38"
    def species = task.ext.species?:"homo_sapiens"
    def merged =  task.ext.merged?:true
    def vep_directory = task.ext.vep_directory?:"/LAB-DATA/GLiCID/projects/BiRD_resources/apps/vep"
    def prefix = task.ext.prefix?:vcf.baseName+".vep"
    def cache_version = task.ext.cache_version?:"113"
"""
mkdir -p TMP


bcftools query -l '${vcf}' |\\
    awk 'END{N=NR;if(N==0) N=1;B=5000.0/N;if(B<1000.0) B=1000;printf("%d\\n",int(B));}' > TMP/jeter.buffer_size


bcftools view "${vcf}" |\\
vep \\
    ${task.cpus>1?"--fork ${task.cpus}":""} \\
    --buffer_size `cat TMP/jeter.buffer_size` \\
    --cache \\
    --dir "${vep_directory}" \\
    --species ${species} \\
    ${merged?"--merged":""} \\
    --assembly ${assembly} \\
    --cache_version "${cache_version}" \\
    --fasta "${params.fasta}" \\
    --offline \\
    --symbol \\
    --format vcf \\
    --force_overwrite \\
    -o TMP/jeter.vcf \\
    --quiet \\
    --vcf \\
    --no_stats

bcftools view --threads ${task.cpus} -O z -o "${prefix}.vcf.gz" TMP/jeter.vcf
bcftools index --threads ${task.cpus} -t -f "${prefix}.vcf.gz"
"""
}
