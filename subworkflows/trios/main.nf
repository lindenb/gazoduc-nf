include {BCTOOLS_MENDELIAN2  } from '../../../modules/bcftools/mendelian2/main.nf'
include {GATK_POSSIBLE_DENOVO} from '../../../modules/gatk/possibledenovo/main.nf'

workflow TRIOS {
    take:
        meta
        fasta
        fai
        dict
        vcf
        pedigree
    main:
        BCTOOLS_MENDELIAN2(fai, pedigree, vcf)
        GATK_POSSIBLE_DENOVO(fasta,fai,dict,BCTOOLS_MENDELIAN.out.vcf,pedigree)
    emit:
        vcf = GATK_POSSIBLE_DENOVO.out.vcf
}





