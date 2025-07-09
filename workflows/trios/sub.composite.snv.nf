include {BCFTOOL_CONCAT                           } from '../../modules/bcftools/concat/main.nf'
include {HET_COMPOSITE                            } from '../../modules/jvarkit/hetcomposite/main.nf'

workflow WORKFLOW_COMPOSITE_SNV {
take:
    meta
    fasta
    fai
    dict
    gff3
    gtf
    triosbams_ch
    pedigree
    vcf
main:
    versions = Channel.empty()

     BCFTOOL_CONCAT(
        vcf
            .map{[it[1],it[2]]}
            .collect()
            .map{[meta,it.flatten()]},
        [[id:"nobed"],[]
        ])
    versions =  versions.mix(BCFTOOL_CONCAT.out.versions)
    vcf = BCFTOOL_CONCAT.out.vcf

    HET_COMPOSITE(
        fasta,fai,dict,
        pedigree,
        vcf
    )
    versions =  versions.mix(HET_COMPOSITE.out.versions)
emit:
    versions
}