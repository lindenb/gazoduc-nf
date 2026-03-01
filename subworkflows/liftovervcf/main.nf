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
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include { DICT_TO_BED                              } from '../../modules/jvarkit/dict2bed'
include { VCF_TO_CONTIGS                           } from '../../subworkflows/bcftools/vcf2contigs'
include { DOWNLOAD_CHAIN                           } from '../../modules/ucsc/download.chain'
include {LIFTOVER_VCF  as LIFT_VCF                 } from '../../modules/gatk/liftovervcf'
include { MERGE_VCFS  as MERGE_VCFS_LIFTED         } from '../../modules/gatk/mergevcfs'
include { MERGE_VCFS  as MERGE_VCFS_FAIL           } from '../../modules/gatk/mergevcfs'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { verify                                   } from '../../modules/utils/functions.nf'

String toUcsc(String build) {
	if(build.equalsIgnoreCase("grch38")) return "hg38";
    if(build.equalsIgnoreCase("grch37")) return "hg19";
	return build;
	}


workflow LIFTOVER_VCF {
    take:
        metadata
        fasta
        fai
        dict
        vcfs //[meta,vcf,tbi]
    main:
        versions = Channel.empty()
        multiqc = Channel.empty()

    	DICT_TO_BED(vcfs.map{meta,vcf,_tbi->[meta,vcf]})
	    versions = versions.mix(DICT_TO_BED.out.versions)
        
        ch1 = DICT_TO_BED.out.bed
            .splitCsv(header:true,sep:'\t')
            .map{meta,row->[meta,row.buildName]}
            .filter{_meta,build->!(isBlank(build) || build==".")}
            .map{meta,build->[meta,toUcsc(build)]}
            .unique()
            .groupTuple()
            .map{meta,builds->{
                verify(builds.size()==1,"multiple builds for ${meta}")
                return [meta,builds[0]];
                }}
            .map{meta,build->[meta,build]}
            .join(vcfs)
            .combine(dict)
            .multiMap{meta1,build,vcf,_tbi,meta2,dict->
                source: [meta1.plus(ucsc_name:build),vcf]
                dest: [meta2,dict]
                }
        DOWNLOAD_CHAIN(ch1.source,ch1.dest)
        versions = versions.mix(DOWNLOAD_CHAIN.out.versions)
        

        VCF_TO_CONTIGS(metadata,vcfs)
        versions = versions.mix(VCF_TO_CONTIGS.out.versions)


        ch2 = DOWNLOAD_CHAIN.out.chain
            .map{meta1,_meta2,chain->[meta1,chain]}
            .combine(VCF_TO_CONTIGS.out.vcf)
            .filter{meta1,_chain,meta2,_vcf,_tbi->meta1.id==meta2.id}
            .multiMap{meta1,liftchain,meta2,vcffile,tbi->
                chain: [meta1,liftchain]
                vcf: [meta2,vcffile,tbi]
                }
        
        LIFT_VCF(
            fasta,
            fai,
           dict,
            ch2.chain,
            ch2.vcf
            )
        versions = versions.mix(LIFT_VCF.out.versions)

        MERGE_VCFS_LIFTED(
            LIFT_VCF.out.vcf
                .map{meta,vcf,tbi->[[id:meta.id],[vcf,tbi]]}
                .groupTuple()
                .map{meta,files->[meta,files.flatten().sort()]}
            )
        versions = versions.mix(MERGE_VCFS_LIFTED.out.versions)

        MERGE_VCFS_FAIL(
            LIFT_VCF.out.fail
                .map{meta,vcf,tbi->[[id:meta.id],[vcf,tbi]]}
                .groupTuple()
                .map{meta,files->[meta,files.flatten().sort()]}
            )
        versions = versions.mix(MERGE_VCFS_FAIL.out.versions)
    emit:
        vcf = MERGE_VCFS_LIFTED.out.vcf
        fail = MERGE_VCFS_FAIL.out.vcf
        versions
        multiqc
}