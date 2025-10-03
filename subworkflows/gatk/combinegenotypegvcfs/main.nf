include {COMBINE_GVCFS     } from '../combinegvcfs'
include {GENOTYPEGVCFS     } from '../../../modules/gatk/genotypegvcfs'

workflow COMBINE_GENOTYPE_GVCFS {
    take:
        meta
        fasta
        fai
        dict
        dbsnp // meta, vcf,tbi
        gvcfs_bed // [meta, gvcfgz, gvcfgz_tbi, bed]
    main:
        versions = Channel.empty()
        COMBINE_GVCFS(meta,fasta,fai,dict,gvcfs_bed)
        versions = versions.mix(COMBINE_GVCFS.out.versions)


        GENOTYPEGVCFS(fasta,fai,dict,dbsnp, COMBINE_GVCFS.out.gvcf)
        versions = versions.mix(GENOTYPEGVCFS.out.versions)
    emit:
        versions
        vcf = GENOTYPEGVCFS.out.vcf

}
