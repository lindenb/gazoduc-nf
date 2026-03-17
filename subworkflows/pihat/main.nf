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
include { isBlank                       } from '../../modules/utils/functions.nf'
include { parseBoolean                  } from '../../modules/utils/functions.nf'
include { VCF_TO_BED                    } from '../../modules/bcftools/vcf2bed'
include { DOWNLOAD_1KG_SAMPLE2POP       } from '../../modules/pihat/download.1kg.pop'
include { PER_CONTIG                    } from '../../modules/pihat/per_contig'
include { PLOT_ASSOC                    } from '../../modules/plink/plot.assoc'
include { PLOT_MDS                      } from '../../modules/plink/plot.mds'
include { PLOT_PIHAT                    } from '../../modules/plink/plot.pihat'
include { AVERAGE_PIHAT                 } from '../../modules/plink/average.pihat'
include { PLINK_ASSOC                   } from '../../modules/plink/assoc'
include { BEDTOOLS_MERGE                } from '../../modules/bedtools/merge'
include { DOWNLOAD_HIGH_LD              } from '../../modules/pihat/high_ld'
include { PLINK_MERGE_BIM_BED_FAM       } from '../../modules/plink/merge'
include { PLINK_GENOME                  } from '../../modules/plink/genome'
include { PLINK_MAKEBED                 } from '../../modules/plink/makebed'
include { PLINK_RECODE_VCF              } from '../../modules/plink/plink2vcf'
include { PLINK_MDS                     } from '../../modules/plink/mds'
include { flatMapByIndex                } from '../../modules/utils/functions.nf'
include { META_TO_PED                   } from '../../subworkflows/pedigree/meta2ped'

String normContig(String s) {
    if(s==null) return "";
    if(!s.matches("(chr)?[0-9]+")) return "";
    if(s.startsWith("chr")) s=s.substring(3);
    return s;
    }

workflow PIHAT {
    take:
        workflow_metadata
        fasta
        fai
        dict
        meta_samples // metas about samples....
        exclude_samples
        exclude_bed
        vcf1kg //[meta,vcf,vcfidx]
        gnomad //[meta,vcf,vcfidx]
        vcfs //[meta,vcf,vcfidx]
    main:


        versions = Channel.empty()
        multiqc = Channel.empty()

        vcfs =  vcfs.map{meta,vcf,tbi->[meta.plus(source:"user"),vcf,tbi]}
        gnomad = gnomad.map{meta,vcf,tbi->[meta.plus(source:"gnomad"),vcf,tbi]}
        vcf1kg = vcf1kg.map{meta,vcf,tbi->[meta.plus(source:"1kg"),vcf,tbi]}
        
        all_kind_of_vcf = vcfs.mix(gnomad).mix(vcf1kg);

        /** check vcf meta.id are unique */
        all_kind_of_vcf
            .map{meta,vcf,tbi->[meta.id,vcf]}
            .groupTuple()
            .filter{id,array->array.size()>1}
            .map{id,array->
                throw new IllegalArgumentException("Duplicate vcf id ${id} for :\n\t${array.join("\n\t")}");
                return -1;
            }

        VCF_TO_BED(all_kind_of_vcf)
        versions = versions.mix(VCF_TO_BED.out.versions)
        
        ch1  = VCF_TO_BED.out.bed.splitCsv(sep:'\t',header:false,elem:1)
            .map{meta,bedrecord->[meta,bedrecord[0]]}
            .combine(all_kind_of_vcf)
            .filter{meta1,_contig,meta2,_vcf,_tbi->meta1.id==meta2.id && meta1.source==meta2.source}
            .map{meta1,contig,_meta2,vcf,tbi->[
                normContig(contig),
                meta1.plus([
                    id : meta1.id+ "."+contig,
                    contig: contig,
                    norm_contig: normContig(contig),
                    ]),
                vcf,
                tbi
                ]}
             .filter{ctg,_meta,_vcf,_tbi->!isBlank(ctg)}



        user_ch = ch1.filter{_ctg,meta,_vcf,_tbi->meta.source=="user"}
        gnomad_ch = ch1.filter{_ctg,meta,_vcf,_tbi->meta.source=="gnomad"}
        k1g_ch = ch1.filter{_ctg,meta,_vcf,_tbi->meta.source=="1kg"}


        user_ch.count()
            .filter{it==0L}
            .map{log.warn("Warning. No chromosome looks like /(chr)?[0-9]+/");}

        dispatch_ch = user_ch.join(k1g_ch,remainder:true)
            .filter{row->row[1]!=null}/* chromosome in 1Kg but not in user vcf */
            .map{row->[
                row[0],/* contig */
                row[1],/* user meta */
                row[2],/* user vcf */
                row[3],/* user tbi */
                (row[4]==null?([:]):row[4]),/* 1kg meta */
                (row[4]==null?([]):row[5]),/* 1kg vcf */
                (row[4]==null?([]):row[6]) /* 1kg tbi */
                ]}.join(gnomad_ch,remainder:true)
                    .filter{row->row[1]!=null}/* chromosome in gnomad but not in user vcf */
                    .map{row->[
                        row[0],/* contig */
                        row[1],/* user meta */
                        row[2],/* user vcf */
                        row[3],/* user tbi */
                        row[4],/* 1kg meta */
                        row[5],/* 1kg vcf */
                        row[6],/* 1kg tbi */
                        (row[7]==null?([:]):row[7]),/* gnomad meta */
                        (row[7]==null?([]):row[8]),/* gnomad vcf */
                        (row[7]==null?([]):row[9]) /* gnomad tbi */
                        ]}.multiMap{_ctg,metau,vcfu,tbiu, metak, vcfk, tbik, metag, vcfg, tbig->
                            user : [metau,vcfu,tbiu]
                            k1g : [metak,vcfk,tbik]
                            gnomad : [metag,vcfg,tbig]
                        }

        /** regions always with high LD */
        DOWNLOAD_HIGH_LD(dict)
        versions = versions.mix(DOWNLOAD_HIGH_LD.out.versions)

        /** merge with user exclude regions */
        BEDTOOLS_MERGE(DOWNLOAD_HIGH_LD.out.bed.mix(exclude_bed)
            .map{_meta,f->f}
            .collect()
            .map{files->[[id:"exclude"],files.sort()]}
            )
        versions = versions.mix(BEDTOOLS_MERGE.out.versions)

	    /* download sample of 1000genome only if vcf1kg defined */
        DOWNLOAD_1KG_SAMPLE2POP(
            k1g_ch
                .count()
                .combine(Channel.of(workflow_metadata))
                .filter{count,_meta-> count>0}
                .map{_count,meta->meta}
                .first()
            )
        versions = versions.mix(DOWNLOAD_1KG_SAMPLE2POP.out.versions)
        
        meta_samples = meta_samples.mix(
            DOWNLOAD_1KG_SAMPLE2POP.out.tsv
                .map{meta,f->f}
                .splitCsv(header:true,sep:'\t')
                .map{meta->[id:meta.Individual_ID, collection:meta.Population, population:meta.Population,  sex:meta.Gender]}
            )

        META_TO_PED(
            workflow_metadata,
            Channel.empty(),
            meta_samples
            )
        versions = versions.mix(META_TO_PED.out.versions)


        PER_CONTIG(
            fasta,
            fai,
            dict,
            exclude_samples.first(),
            BEDTOOLS_MERGE.out.bed,
            dispatch_ch.k1g,
            dispatch_ch.gnomad,
            dispatch_ch.user
            )
        versions = versions.mix(PER_CONTIG.out.versions)
        
     

        /*****************************************************************************/

        PLINK_MERGE_BIM_BED_FAM(
            PER_CONTIG.out.bfile
                .map{_meta,bim,bed,fam->[bim,bed,fam]}
                .flatMap()
                .collect()
                .map{files->[[id:workflow_metadata.id],files.sort()]}
            )
        versions = versions.mix(PLINK_MERGE_BIM_BED_FAM.out.versions)
      
        

        /*****************************************************************************/
	

        PLINK_GENOME(PLINK_MERGE_BIM_BED_FAM.out.bfile)
        versions = versions.mix(PLINK_GENOME.out.versions)
        multiqc = multiqc.mix(PLINK_GENOME.out.related)
        

        bfile_ch = PLINK_MERGE_BIM_BED_FAM.out.bfile
        /**
         * FILTER out related if needed
         */ 
        if(workflow_metadata.filter_out_related==null || parseBoolean(workflow_metadata.filter_out_related)) {
            PLINK_MAKEBED(
                PLINK_GENOME.out.related,
                PLINK_MERGE_BIM_BED_FAM.out.bfile
                )
            versions = versions.mix(PLINK_GENOME.out.versions)
            bfile_ch=  PLINK_MAKEBED.out.bfile
            }

        /**
         *
         * invoke PLINK MDS
         */ 
        PLINK_MDS(
            PLINK_GENOME.out.genome,
            bfile_ch
            )
        versions = versions.mix(PLINK_MDS.out.versions)



        // header of MDS looks like ' FID IID SOL C1 C2 C3 ' . get number of components
        components0_ch = PLINK_MDS.out.mds.splitCsv(header:false,limit:1,sep:' ',strip:true)
            .map{meta,array->array.join(" ")}
            .map{S->S.trim()}
            .map{S->java.util.Arrays.asList(S.split("[ ]+"))}
            .map{array->array.subList(3,array.size())}
            .flatMap()
        
        // all pairs C1,C2 / C1,C3 / C2,C3 etc...
        components1_ch = components0_ch
            .combine(components0_ch)
            .filter{C1,C2 -> C1.compareTo(C2)<0}

        PLOT_MDS(
            META_TO_PED.out.sample2group.first(),
            PLINK_MDS.out.mds
                .combine(components1_ch)
                .map{meta,mds,CX,CY->[meta.plus(Cx:CX,Cy:CY,id:meta.id+".${CX}_${CY}"),mds]}
            )
        versions = versions.mix(PLOT_MDS.out.versions)
        multiqc = multiqc.mix(PLOT_MDS.out.png)

        PLOT_PIHAT(PLINK_GENOME.out.genome)
        versions = versions.mix(PLOT_PIHAT.out.versions)
        multiqc = multiqc.mix(PLOT_PIHAT.out.png)
        multiqc = multiqc.mix(PLOT_PIHAT.out.high_txt)


        AVERAGE_PIHAT(
            PLINK_GENOME.out.genome,
            META_TO_PED.out.sample2group.first()
            )
	    multiqc = multiqc.mix(AVERAGE_PIHAT.out.png)
        multiqc = multiqc.mix(AVERAGE_PIHAT.out.exclude)
        versions = versions.mix(AVERAGE_PIHAT.out.versions)


        PLINK_ASSOC(
            PLINK_MDS.out.mds,
            bfile_ch
            )
        versions = versions.mix(PLINK_ASSOC.out.versions)
        

        PLOT_ASSOC(
            fai,
            PLINK_ASSOC.out.assoc
                .flatMap{row->flatMapByIndex(row,1)}
                .map{meta,assoc->[meta.plus(id:assoc.name),assoc]}
            )
        versions = versions.mix(PLOT_ASSOC.out.versions)
        multiqc = multiqc.mix(PLOT_ASSOC.out.png)
        
    emit:
        versions
        genome = PLINK_GENOME.out.genome
        mds = PLINK_MDS.out.mds
        multiqc
}


