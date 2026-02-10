/*

Copyright (c) 2026 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



*/
include { SPLIT_N_VARIANTS                         } from '../../modules/jvarkit/splitnvariants'
include { WORKFLOW_DENOVO_SNV                      } from './sub.denovo.snv.nf'
include { WORKFLOW_COMPOSITE_SNV                   } from './sub.composite.snv.nf'
include { ANNOT_SNV                                } from  '../../subworkflows/bigannotsnv'
include { flatMapByIndex                           }from '../../modules/utils/functions.nf'

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
    multiqc = Channel.empty()


    

    /** break the VCF into contig */
    SPLIT_N_VARIANTS(
        [[id:"nobed"],[]],
        fai.map{_meta,fai->fai}
            .splitCsv(header:false, sep:'\t')
            .map{row->row[0]}/* extract chromosomes in fai */
            .filter{contig->contig.matches("(chr)?[0-9XY]+")}
            .combine(vcf)
            .map{chr,meta,vcf,tbi->[meta.plus([id:"${meta.id}.${chr}",contig:"${chr}"]),vcf,tbi]}
        )

    
   

    ANNOT_SNV(
        meta,
        fasta,
        fai,
        dict,
        pedigree,
        gtf,
        gff3,
        SPLIT_N_VARIANTS.out.vcf.flatMap(row->flatMapByIndex(row,1))
            .combine(SPLIT_N_VARIANTS.out.tbi.flatMap(row->flatMapByIndex(row,1)))
            .filter{meta1,vcf,meta2,tbi->meta1.id==meta2.id && "${vcf.name}.tbi" == tbi.name}
            .map{meta1,vcf,meta2,tbi->[meta1,vcf,tbi]}
        )

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
        ANNOT_SNV.out.vcf
        )
    versions = versions.mix(WORKFLOW_DENOVO_SNV.out.versions)
    multiqc = multiqc.mix(WORKFLOW_DENOVO_SNV.out.multiqc)
/*
    WORKFLOW_COMPOSITE_SNV(
        meta,
        fasta,
        fai,
        dict,
        gff3,
        gtf,
        triosbams_ch,
        pedigree,
         ANNOT_SNV.out.vcf
        )
     versions = versions.mix(WORKFLOW_COMPOSITE_SNV.out.versions)
     */
emit:
    versions
    multiqc
}
