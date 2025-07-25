include {SPLIT_VCF                                } from '../../subworkflows/jvarkit/splitnvariants/main.nf'
include {WORKFLOW_DENOVO_SNV                      } from './sub.denovo.snv.nf'
include {WORKFLOW_COMPOSITE_SNV                   } from './sub.composite.snv.nf'
include {ANNOT_SNV                                } from  '../../subworkflows/bigannotsnv'

workflow WORKFLOW_SNV {
take:
    meta
    fasta
    fai
    dict
    bed
    pedigree
    gnomad
    gff3
    gtf
    triosbams_ch
    vcf
main:
    versions = Channel.empty()

    /** break the VCF into parts */
    SPLIT_VCF(
        meta,
        fasta,
        fai,
        dict,
        bed,
        vcf
        )
    vcf = SPLIT_VCF.output.vcf


    ANNOT_SNV(
        meta,
        fasta,
        fai,
        dict,
        pedigree,
        gtf,
        gff3,
        vcf
        )
    vcf = ANNOT_SNV.out.vcf
    versions = versions.mix(ANNOT_SNV.out.versions)

    WORKFLOW_DENOVO_SNV(
        meta,
        fasta,
        fai,
        dict,
        gff3,
        gtf,
        triosbams_ch,
        pedigree,
        vcf
        )
    versions = versions.mix(WORKFLOW_DENOVO_SNV.out.versions)
    

    WORKFLOW_COMPOSITE_SNV(
        meta,
        fasta,
        fai,
        dict,
        gff3,
        gtf,
        triosbams_ch,
        pedigree,
        vcf
    )
     versions = versions.mix(WORKFLOW_COMPOSITE_SNV.out.versions)
emit:
    versions
}
