include {isGRCh38} from '../../../modules/utils/k1.nf'
include {VEP as VEP_GRCH38} from './grch38.nf'

workflow VEP {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        versions = Channel.empty()
        if(isGRCh38(fai[1])) {
            VEP_GRCH38(meta,fasta,fai,dict,vcf)
            versions = versions.mix(VEP_GRCH38.out.versions)
            
            vcf = VEP_GRCH38.out.vcf
        } else {
            throw new IllegalArgumentException("VEP: unknown build");
        }
    emit:
        versions
        vcf
}