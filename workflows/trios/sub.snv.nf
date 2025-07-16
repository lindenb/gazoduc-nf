include {SNPEFF as APPLY_SNPEFF                   } from '../../subworkflows/snpeff/main.nf'
include {JVARKIT_VCFGNOMAD                        } from '../../modules/jvarkit/vcfgnomad/main.nf'
include {BCFTOOLS_NORM                            } from '../../modules/bcftools/norm/main.nf'
include {SPLIT_VCF                                } from '../../subworkflows/jvarkit/splitnvariants/main.nf'
include {WORKFLOW_DENOVO_SNV                      } from './sub.denovo.snv.nf'
include {WORKFLOW_COMPOSITE_SNV                   } from './sub.composite.snv.nf'
include {CLINVAR                                  } from '../../subworkflows/annotation/clinvar/main.nf'
include {ALPHAMISSENSE                            } from '../../subworkflows/annotation/alphamissense/main.nf'
include {BHFUCL                                   } from '../../subworkflows/annotation/bhfucl/main.nf'
include {BCFTOOLS_BCSQ                            } from '../../modules/bcftools/bcsq/main.nf'
include {REVEL                                    } from '../../subworkflows/annotation/revel/main.nf'
include {CARDIOPANEL                              } from '../../subworkflows/annotation/cardiopanel/main.nf'
include {PANMASK                                  } from '../../subworkflows/annotation/panmask/main.nf'
include {TISSUES                                  } from '../../subworkflows/jensenlab/tissues/main.nf'
include {DISEASES                                 } from '../../subworkflows/jensenlab/diseases/main.nf'

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


    /** normalize variants */
    BCFTOOLS_NORM(
        fasta,fai,
        [[:],[]],/* no bed */
        vcf
        )
    versions = versions.mix(BCFTOOLS_NORM.out.versions)
    vcf = BCFTOOLS_NORM.out.vcf

    /** filter variants for prediction */
    APPLY_SNPEFF(meta,fasta,fai,dict,vcf)
    versions = versions.mix(APPLY_SNPEFF.out.versions)
    vcf = APPLY_SNPEFF.out.vcf
  
    /** filter variants for gnomad */
    JVARKIT_VCFGNOMAD(gnomad,vcf)
    versions = versions.mix(JVARKIT_VCFGNOMAD.out.versions)
    vcf = JVARKIT_VCFGNOMAD.out.vcf

    BCFTOOLS_BCSQ(fasta,fai,gff3,vcf )
    vcf = BCFTOOLS_BCSQ.out.vcf
    
    CLINVAR(meta,fasta,fai,dict,[[id:"nobed"],[]],vcf)
    vcf = CLINVAR.out.vcf

    ALPHAMISSENSE(meta,fasta,fai,dict,[[id:"nobed"],[]],vcf)
    vcf = ALPHAMISSENSE.out.vcf

    BHFUCL(meta,fasta,fai,dict,gtf,vcf)
    vcf = BHFUCL.out.vcf
 
    REVEL(meta,fasta,fai,dict,vcf)
    vcf = REVEL.out.vcf

    CARDIOPANEL(meta,fasta,fai,dict,gtf,vcf)
    vcf = CARDIOPANEL.out.vcf

    PANMASK(meta, fasta, fai, dict, vcf)
    versions = versions.mix(PANMASK.out.versions)
    vcf = PANMASK.out.vcf

    TISSUES(meta, fasta, fai, dict, gtf, vcf)
	versions = versions.mix(TISSUES.out.versions)
	vcf = TISSUES.out.vcf


    DISEASES(meta, fasta, fai, dict, gtf, vcf)
	versions = versions.mix(DISEASES.out.versions)
	vcf = DISEASES.out.vcf


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
