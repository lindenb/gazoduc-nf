workflow TRIOS {
    take:
        meta
        fasta
        fai
        dict
        vcf
        pedigree
    main:
      BCTOOLS_MENDELIAN(fai,pedigree,vcf)

}


process BCTOOLS_MENDELIAN {
    tag "${contig}"
    input:
        tuple val(meta1),path(fai)
        tuple val(meta),val(contig),path(vcfidx)
        tuple val(meta2),path(pedigree)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def prefix = meta.id+".md."+contig
    """
    bcftools +mendelian2 -O z -o "${prefix}.vcf.gz" --write-index --rules "\${URL}" --regions "${contig}" -m a -P '${pedigree}' "${vcf}"
    """
    }


process GATK_POSSIBLE_DENOVO {
    tag "${contig}"
    input:
        tuple val(meta),path(vcf),path(vcfidx)
        tuple val(meta2),path(pedigree)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def prefix = meta.id+".md."+contig
    """
    gatk --java-options "-Xmx${task.memory.giga}g" VariantAnnotator \
        --annotation PossibleDeNovo \
        --pedigree  '${pedigree}' \\
        -V "${vcf}" \\
        -O "${prefix}.vcf.gz"
    """
    }


