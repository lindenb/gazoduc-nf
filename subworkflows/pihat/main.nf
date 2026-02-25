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
include { VCF_TO_BED                    } from '../../modules/bcftools/vcf2bed'
include { DOWNLOAD_1KG_SAMPLE2POP       } from '../../modules/pihat/download.1kg.pop'
include { PER_CONTIG                    } from '../../modules/pihat/per_contig'
include { PLOT_ASSOC                    } from '../../modules/pihat/plot.assoc'
include { PLOT_MDS                      } from '../../modules/pihat/plot.mds'
include { PLOT_PIHAT                    } from '../../modules/pihat/plot.pihat'
include { AVERAGE_PIHAT                 } from '../../modules/pihat/average.pihat'
include { PLINK_ASSOC                   } from '../../modules/pihat/plink.assoc'
include { PLINK_GENOME                  } from '../../modules/pihat/plink.genome'
include { BEDTOOLS_MERGE                } from '../../modules/bedtools/merge'
include { DOWNLOAD_HIGH_LD              } from '../../modules/pihat/high_ld'

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
        sample2pop
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

        VCF_TO_BED(all_kind_of_vcf)
        versions = versions.mix(VCF_TO_BED.out.versions)
        
        ch1  = VCF_TO_BED.out.bed.splitCsv(sep:'\t',header:false,elem:1)
            .map{meta,bedrecord->[meta,bedrecord[0]]}
            .combine(all_kind_of_vcf)
            .filter{meta1,contig,meta2,_vcf,_tbi->meta1.id==meta2.id && meta1.source==meta2.source}
            .map{meta1,contig,meta2,vcf,tbi->[
                normContig(contig),
                meta1.plus([
                    id : meta1.id+ "."+contig,
                    contig: contig,
                    norm_contig: normContig(contig),
                    ]),
                vcf,
                idx
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
        BEDTOOLS_MERGE(DOWNLOAD_HIGH_LD.out.mix(exclude_bed)
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

        
        PLINK_GENOME(
            fasta,
            fai,
            dict,
            PER_CONTIG.out.bfile
                .collect()
                .map{[[id:workflow_metadata.id],it]}
            )
        versions = versions.mix(PLINK_GENOME.out.versions)


        MERGE_SAMPLE2POP(
            DOWNLOAD_1KG_SAMPLE2POP.out.tsv,
            sample2pop,
            PLINK_GENOME.out.merged_plink
            )
        versions = versions.mix(MERGE_SAMPLE2POP.out.versions)


        PLINK_ASSOC(fasta, fai, dict, PLINK_GENOME.out.merged_plink , PLINK_GENOME.out.mds )
        versions = versions.mix(PLINK_ASSOC.out.versions)


        PLOT_ASSOC(
            fasta,
            fai,
            dict,
            PLINK_ASSOC.out.assoc
                .flatMap{it[1]}
                .filter{it.name.matches(".*C[123].qassoc")}
                .map{[[id:"pihat"],it]}
            )
        versions = versions.mix(PLOT_ASSOC.out.versions)

        PLOT_PIHAT(PLINK_GENOME.out.genome)
        versions = versions.mix(PLOT_PIHAT.out.versions)
	    multiqc = multiqc.mix(PLOT_PIHAT.out.pict.filter{_meta,img->img.name.endsWith(".png")})

        components = Channel.of(["C1","C2"],["C1","C3"],["C2","C3"])
        formats = Channel.of("pdf","png")
       
        PLOT_MDS(PLINK_GENOME.out.mds.combine(components).combine(formats))
        versions = versions.mix(PLOT_MDS.out.versions)
       
	multiqc = multiqc.mix(PLOT_MDS.out.plot.filter{_meta,img->img.name.endsWith(".png")})

        AVERAGE_PIHAT(PLINK_GENOME.out.genome, MERGE_SAMPLE2POP.out.sample2pop)
	multiqc = multiqc.mix(AVERAGE_PIHAT.out.png)
        versions = versions.mix(AVERAGE_PIHAT.out.versions)
    emit:
        versions
        genome = PLINK_GENOME.out.genome
        mds = PLINK_GENOME.out.mds
	multiqc
}


process MERGE_SAMPLE2POP {
tag "${meta1.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
input:
    tuple val(meta1),path(opt_tsv1)
    tuple val(meta2),path(opt_tsv2)
    tuple val(meta),path(merged_plink)
output:
    tuple val(meta1),path("*.tsv"),emit:sample2pop
    path("versions.yml"),emit:versions
script:
    def fam = merged_plink.find{it.name.endsWith(".fam")}
    def prefix = task.ext.prefix?:"merged_sample2pop"
    def other_name = task.ext.other?:"OTHER"
"""
touch jeter.tsv

if ${opt_tsv1?true:false}
then
    cat ${opt_tsv1} >> jeter.tsv
fi

if ${opt_tsv2?true:false}
then
    cat ${opt_tsv2} >> jeter.tsv
fi

cut -f1,2 jeter.tsv |\\
    sort -T . -t '\t' -k1,1 --unique > jeter2.tsv

mv jeter2.tsv jeter.tsv

# all uniq names in 1kg and sample2pop
cut -f1 jeter.tsv | sort -T TMP | uniq >  jeter2.tsv
# sample names in plink.fam
awk '{print \$2}' "${fam}" | sort -T TMP | uniq >  jeter3.tsv
# extract sample without name in jeter.tsv
comm -13 jeter2.tsv jeter3.tsv |\\
    awk '{printf("%s\t${other_name}\\n",\$1);}' >> jeter.tsv

rm jeter2.tsv jeter3.tsv
mv jeter.tsv '${prefix}.tsv'

cat << EOF > versions.yml
${task.process}:
    sort: todo
EOF
"""
}







