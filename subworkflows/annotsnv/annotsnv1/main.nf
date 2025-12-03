/*

Copyright (c) 2025 Pierre Lindenbaum

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
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include { isBlank             } from '../../../modules/utils/functions'
include { parseBoolean        } from '../../../modules/utils/functions'
include { SNPEFF              } from '../../../subworkflows/snpeff'
include { JVARKIT_VCFGNOMAD   } from '../../../modules/jvarkit/vcfgnomad'
include { JVARKIT_VCFFILTERSO } from '../../../modules/jvarkit/vcffilterso'
include { BCFTOOLS_NORM       } from '../../../modules/bcftools/norm'

workflow ANNOT_SNV {
take:
    metadata
    fasta
    fai
    dict
    snpeff_directory
    gnomad_vcf
    gtf
    gff3
    pedigree
    vcf
main:
    versions = Channel.empty()
    multiqc = Channel.empty()

    if(metadata.with_snpeff==null) {
        log.warn("ANNOT_SNV  : with_snpeff undefined")
        }  
    if(metadata.with_gnomad==null) {
        log.warn("ANNOT_SNV  : with_gnomad undefined")
        } 
    if(metadata.bcftools_norm==null) {
        log.warn("ANNOT_SNV  : with_bcftools_norm undefined")
        }
     if(metadata.with_filterso==null) {
        log.warn("ANNOT_SNV  : with_filterso undefined")
        } 
    
    if(metadata.bcftools_norm!=null && parseBoolean(metadata.bcftools_norm)) {
        BCFTOOLS_NORM(fasta,fai,vcf)
        versions = versions.mix(BCFTOOLS_NORM.out.versions)
        vcf = BCFTOOLS_NORM.out.vcf
        }
    
    if(metadata.with_snpeff==null || parseBoolean(metadata.with_snpeff)) {
        SNPEFF(metadata,snpeff_directory,vcf)
        versions = versions.mix(SNPEFF.out.versions)
        multiqc = versions.mix(SNPEFF.out.multiqc)
        vcf = SNPEFF.out.vcf
        }
    
     if(metadata.with_filterso==null || parseBoolean(metadata.with_filterso)) {
        JVARKIT_VCFFILTERSO(vcf)
        versions = versions.mix(JVARKIT_VCFFILTERSO.out.versions)
        vcf = JVARKIT_VCFFILTERSO.out.vcf
        }

    if(metadata.with_gnomad==null || parseBoolean(metadata.with_gnomad)) {
        JVARKIT_VCFGNOMAD(
            gnomad_vcf,
            vcf
            )
        versions = versions.mix(JVARKIT_VCFGNOMAD.out.versions)
        vcf = JVARKIT_VCFGNOMAD.out.vcf
        }
emit:
    vcf
    versions
    multiqc
}
