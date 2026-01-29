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
nextflow.enable.dsl=2

include { validateParameters                 } from 'plugin/nf-schema'
include { paramsHelp                         } from 'plugin/nf-schema'
include { paramsSummaryLog                   } from 'plugin/nf-schema'
include { samplesheetToList                  } from 'plugin/nf-schema'
include {assertKeyExistsAndNotEmpty          } from '../../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty            } from '../../../modules/utils/functions.nf'
include {assertKeyMatchRegex                 } from '../../../modules/utils/functions.nf'
include {parseBoolean                        } from '../../../modules/utils/functions.nf'
include {isBlank                             } from '../../../modules/utils/functions.nf'
include {VCF_STATS                           } from '../../../subworkflows/vcfstats'
include {BEDTOOLS_MAKEWINDOWS                } from '../../../modules/bedtools/makewindows'
include {BED_CLUSTER                         } from '../../../modules/jvarkit/bedcluster'
include {GATK_BAM2VCF                        } from '../../../subworkflows/gatk/bam2vcf'
include {GATK_BAM2VCF_DNC                    } from '../../../subworkflows/gatk/bam2vcf.dnc'
include {runOnComplete; dumpParams           } from '../../../modules/utils/functions.nf'
include {BCFTOOLS_GUESS_PLOIDY               } from '../../../modules/bcftools/guess_ploidy'
include {BCFTOOLS_STATS                      } from '../../../modules/bcftools/stats'
include {HAPLOTYPECALLER                     } from '../../../subworkflows/gatk/haplotypecaller'
include {HAPLOTYPECALLER_DIRECT              } from '../../../subworkflows/gatk/haplotypecaller.direct'
include {SAMTOOLS_SAMPLES                    } from '../../../modules/samtools/samples'
include {META_TO_PED                         } from '../../../subworkflows/pedigree/meta2ped'
include {PREPARE_ONE_REFERENCE               } from '../../../subworkflows/samtools/prepare.one.ref'
include {MULTIQC                             } from '../../../subworkflows/multiqc'
include {PREPARE_USER_BED                    } from '../../../subworkflows/bedtools/prepare.user.bed'
include {VCF_INPUT                           } from '../../../subworkflows/nf/vcf_input'
include {READ_SAMPLESHEET                    } from '../../../subworkflows/nf/read_samplesheet'
include {META_TO_BAMS                        } from '../../../subworkflows/samtools/meta2bams2'
include { GTF_INPUT                          } from '../../../subworkflows/nf/gtf_input'
include { GTF_TO_EXOME                       } from '../../../modules/gtf/gtf2exome1'
include {JVARKIT_VCFFILTERJDK                } from '../../../modules/jvarkit/vcffilterjdk'
include {ANNOT_SNV                           } from '../../../subworkflows/annotsnv/annotsnv1'
include {BCFTOOLS_SORT                       } from '../../../modules/bcftools/sort'
include {BCFTOOLS_INDEX                      } from '../../../modules/bcftools/index'
include {BCFTOOLS_CONCAT as CONCAT_ANNOT     } from '../../../modules/bcftools/concat3'
include { JVARKIT_VCF_TO_TABLE as VCF2TXT    } from '../../../modules/jvarkit/vcf2table'
include { JVARKIT_VCF_TO_TABLE as VCF2HTML   } from '../../../modules/jvarkit/vcf2table'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline ")
   exit 0
}


workflow {

    if(!workflow.stubRun) {
		validateParameters()
		}
	if( params.help ) {
		log.info(paramsHelp())
            exit 0
        }  else {
        // Print summary of supplied parameters
        log.info paramsSummaryLog(workflow)
        }

    /* no fastq samplesheet */
	if(params.samplesheet==null) {
		log.error("--samplesheet undefined")
		exit -1
		}
	/* no fastq samplesheet */
	if(params.fasta==null) {
		log.error("--fasta undefined")
		exit -1
		}
	bams_ch = Channel.empty()
	versions = Channel.empty()
	multiqc = Channel.empty()

	workflow_metadata = [
		id: "haplotypecaller"
		]



    PREPARE_ONE_REFERENCE(
        workflow_metadata,
        Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
        )
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)




    if(params.dbsnp==null)
        {
        dbsnp =     Channel.of([workflow_metadata, [] ,[] ])
        }
    else
        {
        VCF_INPUT([
            arg_name       : "dbsnp",
            path          : params.dbsnp,
            require_index : true,
            unique        : true
            ])
        versions = versions.mix(VCF_INPUT.out.versions)
        dbsnp = VCF_INPUT.out.vcf
        }

    def gtf     =  Channel.of([workflow_metadata,[],[]]).first()
  
    if(params.gtf==null) {
      /***************************************************
		*
		*  DOWNLOAD GTF
		*
		*/
		GTF_INPUT(
				workflow_metadata.plus([
					arg_name: "gtf",
					require_index: true,
					download: true,
					path: params.gtf
					]),
				PREPARE_ONE_REFERENCE.out.dict
				)
		versions = versions.mix(GTF_INPUT.out.versions)
        gtf = GTF_INPUT.out.gtf
        }


	READ_SAMPLESHEET(
        workflow_metadata.plus(arg_name:"samplesheet"),
        params.samplesheet
        )
     versions = versions.mix(READ_SAMPLESHEET.out.versions)


    META_TO_BAMS(
        workflow_metadata,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        READ_SAMPLESHEET.out.samplesheet
        )
    versions = versions.mix(META_TO_BAMS.out.versions)
  


    bams_ch = META_TO_BAMS.out.bams.map{meta,bam,bai,fasta,fai,dict->[meta,bam,bai]}

    /* check no duplicate samples */
    bams_ch.map{meta,bam,bai->meta.id}
        .unique()
        .count()
        .combine(bams_ch.map{meta,bam,bai->meta.id}.count())
        .filter{c1,c2->c1!=c2}
        .map{
            throw new IllegalArgumentException("Check the samplesheet. There is a duplicate sample name");
            }
    
    /** build pedigree from meta data */
    META_TO_PED(workflow_metadata, bams_ch.map{it[0]})
    versions = versions.mix(META_TO_PED.out.versions)


    if(params.bed==null) {
       bed = PREPARE_ONE_REFERENCE.out.scatter_bed
        }
    else if(params.bed=="exome") {
		GTF_TO_EXOME(
			PREPARE_ONE_REFERENCE.out.fai,
			gtf.map{meta,gtf_file,_tbi->[meta,gtf_file]}
			)
		versions = versions.mix(GTF_TO_EXOME.out.versions)
		bed  = GTF_TO_EXOME.out.bed
        }
    else {
		bed = Channel.of([[id:file(params.bed).baseName],file(params.bed)])
        
        PREPARE_USER_BED(
            workflow_metadata,
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            PREPARE_ONE_REFERENCE.out.scatter_bed,
            bed
            )
        versions = versions.mix(PREPARE_USER_BED.out.versions)
        multiqc = multiqc.mix(PREPARE_USER_BED.out.multiqc)
        bed = PREPARE_USER_BED.out.bed.first()
        }

   /* cut the bed/genome into parts for SV calling per region */
   BEDTOOLS_MAKEWINDOWS(bed)
   versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

   /* if it's an exome , group the small genome together in BED */
   BED_CLUSTER(
        PREPARE_ONE_REFERENCE.out.dict,
        BEDTOOLS_MAKEWINDOWS.out.bed
        )
   versions = versions.mix(BED_CLUSTER.out.versions)
   beds_ch = BED_CLUSTER.out.bed
        .map{_meta,beds->beds}
        .map{beds->beds instanceof List?beds:[beds]}
        .flatMap()
        .map{bed->[[id:bed.baseName],bed]}


    vcf_ch = Channel.empty()

    if(params.method.equalsIgnoreCase("gvcf")) {
        HAPLOTYPECALLER(
            [id:"hapcaller",gvcf_merge_method:params.gvcf_merge_method],
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            META_TO_BAMS.out.all_references,
            dbsnp,
            META_TO_PED.out.pedigree_gatk,
            beds_ch,
            bams_ch
            )
        versions = versions.mix(HAPLOTYPECALLER.out.versions)
        vcf_ch = HAPLOTYPECALLER.out.vcf
        }
    else if(params.method.equalsIgnoreCase("bam2vcf")) {
        if(params.divide_and_conquer==true) {
            GATK_BAM2VCF_DNC(
                workflow_metadata,
                PREPARE_ONE_REFERENCE.out.fasta,
                PREPARE_ONE_REFERENCE.out.fai,
                PREPARE_ONE_REFERENCE.out.dict,
                dbsnp,
                META_TO_PED.out.pedigree_gatk,
                META_TO_BAMS.out.all_references, //[meta, [ref files fa fai dict...]] all known reference
                beds_ch, // [meta,bed]
                bams_ch // [meta,bam,bai]
            )
        } else {
            GATK_BAM2VCF(
                workflow_metadata,
                PREPARE_ONE_REFERENCE.out.fasta,
                PREPARE_ONE_REFERENCE.out.fai,
                PREPARE_ONE_REFERENCE.out.dict,
                dbsnp,
                META_TO_PED.out.pedigree_gatk,
                META_TO_BAMS.out.all_references, //[meta, [ref files fa fai dict...]] all known reference
                beds_ch, // [meta,bed]
                bams_ch // [meta,bam,bai]
                )
            versions = versions.mix(GATK_BAM2VCF.out.versions)
            vcf_ch = GATK_BAM2VCF.out.vcf
            }
        }
    else if(params.method.equalsIgnoreCase("direct")) {
        HAPLOTYPECALLER_DIRECT(
            workflow_metadata,
             PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            dbsnp,
            META_TO_PED.out.pedigree_gatk,
            beds_ch,
            bams_ch
            )
        versions = versions.mix(HAPLOTYPECALLER_DIRECT.out.versions)
        vcf_ch = HAPLOTYPECALLER_DIRECT.out.vcf
        }
    else
        {
        throw new IllegalArgumentException("undefined params.method=${params.method}")
        }
 
 /***************************************************
   *
   * BCFTOOLS_GUESS_PLOIDY
   *
   */
if( parseBoolean(params.with_guess_ploidy)) {
        BCFTOOLS_GUESS_PLOIDY(
                PREPARE_ONE_REFERENCE.out.fasta,
                PREPARE_ONE_REFERENCE.out.fai,
                vcf_ch
                )
        versions = versions.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)
        multiqc = multiqc.mix(BCFTOOLS_GUESS_PLOIDY.out.output)
        }
    /***************************************************
    *
    * JVARKIT FILTERJDK
    *
    */
    if(params.jvarkit_filter!=null) {
        JVARKIT_VCFFILTERJDK(
            jvarkit_filter,
            [[id:"no_ped"],[]],
            vcf_ch
		    )
	    versions = versions.mix(JVARKIT_VCFFILTERJDK.out.versions)
	    vcfs = JVARKIT_VCFFILTERJDK.out.vcf
        }


    if( parseBoolean(params.with_annotation)) {
        if(params.snpeff_database_directory==null) {
            snpeff_database_directory = [[id:"nosnpeffdir"],[]]
            }
        else
            {
            snpeff_database_directory = PREPARE_ONE_REFERENCE.out.fai
                .map{meta,_fai->[meta,file(params.snpeff_database_directory)]}
                
            }
  
        if(params.gnomad==null) {
            gnomad_vcf = [[id:"gnomad"],[],[]]
            }
        else
            {
            gnomad_vcf = [[id:"gnomad"],file(params.gnomad),file(params.gnomad+".tbi")]
            }

        ANNOT_SNV(
            workflow_metadata.plus(
                with_snpeff: (params.snpeff_database_directory!=null),
                with_bcftools_norm : parseBoolean(params.with_bcftools_norm),
                with_filterso : !isBlank(params.soacn),
                with_gnomad : (params.gnomad!=null)
                ),
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            snpeff_database_directory,
            gnomad_vcf,
            [[id:"nogtf"],[]],
            [[id:"nogff3"],[]],
            [[id:"noped"],[]],
            vcf_ch.map{meta,vcf,idx->[meta,vcf]}
            )
        versions = versions.mix(ANNOT_SNV.out.versions)
        multiqc = multiqc.mix(ANNOT_SNV.out.multiqc)
	    vcfs = ANNOT_SNV.out.vcf

        BCFTOOLS_SORT(vcfs)
        versions = versions.mix(BCFTOOLS_SORT.out.versions)
        BCFTOOLS_INDEX(BCFTOOLS_SORT.out.vcf)
        versions = versions.mix(BCFTOOLS_INDEX.out.versions)

        CONCAT_ANNOT(BCFTOOLS_INDEX.out.vcf
                .map{meta,vcf,tbi->[meta.plus(id:"${meta.id}.annot"),[vcf,tbi]]}
                .groupTuple()
                .map{meta,files->[meta,files.flatten().sort()]}
            )
        versions = versions.mix(CONCAT_ANNOT.out.versions)
        vcfs = CONCAT_ANNOT.out.vcf
        }

  /***************************************************
   *
   * VCF TO TABLE
   *
   */
   if( parseBoolean(params.with_vcf2table) ) {
        VCF2TXT(
            [[id:"no_ped"],[]],
            vcfs.map{meta,vcf,_tbi->[meta,vcf]}
            )
        versions = versions.mix(VCF2TXT.out.versions)

        VCF2HTML(
            [[id:"no_ped"],[]],
            vcfs.map{meta,vcf,_tbi->[meta,vcf]}
            )
        versions = versions.mix(VCF2HTML.out.versions)
        }

  /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
   if(parseBoolean(params.with_qc)) {
    BCFTOOLS_STATS(
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        bed,
        gtf.map{[it[0],it[1]]}.first(),//meta,gtf
        [[:],[]],//samples,
        vcf_ch.map{[it[0],[it[1],it[2]]]}
        )
    versions = versions.mix(BCFTOOLS_STATS.out.versions)
    multiqc = multiqc.mix(BCFTOOLS_STATS.out.stats)
    }

 if(parseBoolean(params.with_multiqc)) {
	MULTIQC(
		workflow_metadata.plus("id":"hapcaller"),
		META_TO_PED.out.sample2collection,
		versions,
		[[id:"no_mqc_config"],[]],
		multiqc
		)
    }

}

runOnComplete(workflow)
