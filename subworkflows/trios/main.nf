include {BCTOOLS_MENDELIAN2                 } from '../../modules/bcftools/mendelian2/main.nf'
include {GATK_CALCULATE_GENOTYPE_POSTERIORS } from '../../modules/gatk/calculategenotypeposteriors/main.nf'
include {GATK_POSSIBLE_DENOVO               } from '../../modules/gatk/possibledenovo/main.nf'

workflow TRIOS {
    take:
        meta
        fasta
        fai
        dict
        pedigree
        vcf
    main:
        versions = Channel.empty()
        gatk_supporting = [[id:"empty"], [],[]]
        BCTOOLS_MENDELIAN2(fai, pedigree, vcf)
         versions = versions.mix(BCTOOLS_MENDELIAN2.out.versions)

        GATK_CALCULATE_GENOTYPE_POSTERIORS(fasta,fai,dict,pedigree,gatk_supporting,BCTOOLS_MENDELIAN2.out.vcf)
        versions = versions.mix(GATK_CALCULATE_GENOTYPE_POSTERIORS.out.versions)

        GATK_POSSIBLE_DENOVO(fasta,fai,dict,pedigree,GATK_CALCULATE_GENOTYPE_POSTERIORS.out.vcf)
        versions = versions.mix(GATK_POSSIBLE_DENOVO.out.versions)
    emit:
        vcf = GATK_POSSIBLE_DENOVO.out.vcf
        versions
}




