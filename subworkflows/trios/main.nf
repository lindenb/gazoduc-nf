include {BCTOOLS_MENDELIAN2  } from '../../modules/bcftools/mendelian2/main.nf'
include {GATK_POSSIBLE_DENOVO} from '../../modules/gatk/possibledenovo/main.nf'

workflow TRIOS {
    take:
        meta
        fasta
        fai
        dict
        pedigree
        vcf
    main:
        BCTOOLS_MENDELIAN2(fai, pedigree, vcf)
        GATK_POSSIBLE_DENOVO(fasta,fai,dict,pedigree,BCTOOLS_MENDELIAN2.out.vcf)
    emit:
        vcf = GATK_POSSIBLE_DENOVO.out.vcf
}




