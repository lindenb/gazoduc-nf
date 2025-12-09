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
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



*/

include { validateParameters                      } from 'plugin/nf-schema'
include { paramsHelp                              } from 'plugin/nf-schema'
include { paramsSummaryLog                        } from 'plugin/nf-schema'
include { samplesheetToList                       } from 'plugin/nf-schema'

include {runOnComplete                            } from '../../modules/utils/functions.nf'
include {parseBoolean                             } from '../../modules/utils/functions.nf'
include { isBlank                                 } from '../../modules/utils/functions.nf'
include {FASTP                                    } from '../../subworkflows/fastp'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions'
include {BWA_INDEX                                } from '../../modules/bwa/index'
include { DRAGMAP_INDEX                           } from '../../modules/dragmap/index'
include {PB_FQ2BAM                                } from '../../modules/parabricks/fq2bam'
include {PB_HAPLOTYPECALLER                       } from '../../subworkflows/parabricks/haplotypecaller'
include {PB_DEEPVARIANT                           } from '../../subworkflows/parabricks/deepvariant'
include {PB_BAM2FQ                                } from '../../modules/parabricks/bam2fq'
include {SCATTER_TO_BED                           } from '../../subworkflows/gatk/scatterintervals2bed'
include {BAM_QC                                   } from '../../subworkflows/bamqc'
include {SOMALIER_BAMS                            } from '../../subworkflows/somalier/bams'
include {MANTA                                    } from '../../subworkflows/manta'
include {INDEXCOV                                 } from '../../subworkflows/indexcov'
include {DELLY                                    } from '../../subworkflows/delly2'
include {CNVNATOR                                 } from '../../subworkflows/cnvnator'
include {BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include {BED_CLUSTER as BED_CLUSTER1              } from '../../modules/jvarkit/bedcluster'
include {BED_CLUSTER as BED_CLUSTER2              } from '../../modules/jvarkit/bedcluster'
include {BCFTOOLS_GUESS_PLOIDY                    } from '../../modules/bcftools/guess_ploidy'
include {BCFTOOLS_STATS                           } from '../../modules/bcftools/stats'
include {GRAPHTYPER                               } from '../../subworkflows/graphtyper/genotype'
include {BCFTOOLS_CALL                            } from '../../subworkflows/bcftools/call'
include {FREEBAYES_CALL                           } from '../../subworkflows/freebayes/call'
include {HAPLOTYPECALLER                          } from '../../subworkflows/gatk/haplotypecaller'
include {SAMPLESHEET_TO_FASTQ                     } from '../../subworkflows/samplesheet2fastq'
include {PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include {META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { META_TO_PED                             } from '../../subworkflows/pedigree/meta2ped'
include { PREPARE_USER_BED                        } from '../../subworkflows/bedtools/prepare.user.bed'
include { MAP_BWA                                 } from '../../subworkflows/bwa/map.fastqs'
include { MAP_DRAGMAP                             } from '../../subworkflows/dragmap'
include { NOT_EMPTY_VCF                           } from '../../modules/bcftools/notempty'
include { READ_SAMPLESHEET                        } from '../../subworkflows/nf/read_samplesheet'
include { SMOOVE                                  } from '../../subworkflows/smoove' 
include { VCF_INPUT as DBSNP_VCF_INPUT            } from '../../subworkflows/nf/vcf_input' 
include { VCF_INPUT as KNOWN_INDELS_VCF_INPUT     } from '../../subworkflows/nf/vcf_input'
include { VCF_INPUT as SOMALIER_SITES_VCF_INPUT   } from '../../subworkflows/nf/vcf_input' 
include { GTF_INPUT                               } from '../../subworkflows/nf/gtf_input' 
include { ENCODE_BLACKLIST                        } from '../../modules/encode/blacklist' 
include { BEDTOOLS_SUBTRACT                       } from '../../modules/bedtools/subtract' 
include {UNMAPPED                                 } from '../../subworkflows/unmapped'



workflow {

validateParameters()


if( params.help ) {
    log.info(paramsHelp())
    exit 0
}  else {
  // Print summary of supplied parameters
  log.info paramsSummaryLog(workflow)
}



  if(params.fasta==null) {
      throw new IllegalArgumentException("undefined --fasta");
      }
	
  if(!params.hts_type.matches("WGS|WES")) {
    throw new IllegalArgumentException("hts_type should batch WES|WGS ${params.hts_type}")
    }
  
  def metadata= [
      id: "parabricks"
      ]
  

  def prevent_multi_gpu_ch = Channel.of(1L)



 
  def is_exome = params.hts_type == "WES";
  def is_wgs = params.hts_type == "WGS";

  versions = Channel.empty()
  multiqc  = Channel.empty()

  /* no fastq samplesheet */
  if(params.samplesheet==null) {
    samplesheet0_ch = Channel.empty()
    }
  else
    {
    /* no fastq samplesheet */
    samplesheet0_ch = READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        ).samplesheet
    }
  
  if(!params.mapper.matches("pb|parabricks") && parseBoolean(params.with_merge_fastqs) && parseBoolean(params.with_seqkit_split2)) {
    log.warn("both using --with_merge_fastqs and --with_seqkit_split2 ? ok, why not ?")
    }

  SAMPLESHEET_TO_FASTQ(
    metadata.plus([
      bam2fastq_method: params.bam2fastq_method,
      merge_fastqs : params.mapper.matches("pb|parabricks") || parseBoolean(params.with_merge_fastqs)
      ]),
    samplesheet0_ch
    )
  versions =  SAMPLESHEET_TO_FASTQ.out.versions


  single_end = SAMPLESHEET_TO_FASTQ.out.single_end
    .map{[it[0],[it[1]]]}
    
    
  paired_end = SAMPLESHEET_TO_FASTQ.out.paired_end
    .map{[it[0],[it[1],it[2]]]}

  if(parseBoolean(params.with_fastp)) {
    FASTP(
      metadata.plus([
        fastqc_before : params.with_fastqc,
        fastqc_after : params.with_fastqc,
        fastp_disabled : !(params.with_fastp)
        ]),
      single_end.mix(paired_end) 
      )


    multiqc = multiqc.mix(FASTP.out.multiqc)
    versions = versions.mix(FASTP.out.versions)
    paired_end = FASTP.out.paired_end
    single_end = FASTP.out.single_end
    }
   else
	{
	paired_end = paired_end.map{meta,fqs->[meta,fqs[0],fqs[1]]}
	}

  /***************************************************
   *
   *  DBSNP VCF
   *
   */
  if(params.dbsnp!=null) {
    DBSNP_VCF_INPUT(
      metadata.plus([
        arg_name : "dbsnp",
        path: params.dbsnp,
        require_index : true,
        required : true,
        unique : true
        ])
      )
    versions = versions.mix(DBSNP_VCF_INPUT.out.versions)
    dbsnp_ch = DBSNP_VCF_INPUT.out.vcf.first()
    }
  else
    {
    dbsnp_ch = Channel.of([[id:"nodbsnp"],[],[]]).first()
    }

  /***************************************************
   *
   *  KNWON INDELS
   *
   */
if(params.known_indels_vcf!=null) {
    KNOWN_INDELS_VCF_INPUT(
      metadata.plus([
        arg_name : "known_indels_vcf",
        path: params.known_indels_vcf,
        require_index : true,
        required : true,
        unique : true
        ])
      )
    versions = versions.mix(KNOWN_INDELS_VCF_INPUT.out.versions)
    known_indels_ch = KNOWN_INDELS_VCF_INPUT.out.vcf.first()
    }
  else
    {
    if(parseBoolean(params.with_bqsr)) {
      log.warn("--with_bqsr==true but --known_indels_vcf is not declared.");
      }

    known_indels_ch = Channel.of([[id:"known_indels"],[],[]]).first()
    }


  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


  /***************************************************
   *
   *  GTF
   *
   */
  GTF_INPUT(
    metadata.plus(
      arg_name : "gtf",
      path : params.gtf,
      require_index  : true,
      download: true
      ),
    PREPARE_ONE_REFERENCE.out.dict
    )
  gtf_ch = GTF_INPUT.out.gtf.ifEmpty([[id:"no_gtf"],[],[]])
  

  /***************************************************
   *
   *  Index FASTA For BWA
   *
   */
  if(params.mapper.matches("(pb|parabricks|bwa)")) {
    DragmapDir = Channel.empty()

  	if(params.bwa_index_directory!=null) {
      BWADir = [[id:"bwaindex"],file(params.bwa_index_directory)];
      }
    else
      {
      BWA_INDEX(PREPARE_ONE_REFERENCE.out.fasta)
      versions = versions.mix(BWA_INDEX.out.versions)
      BWADir = BWA_INDEX.out.bwa_index
      }
    }
  /***************************************************
   *
   *  INDEX FASTA FOR DRAGMAP
   *
   */
  else if(params.mapper.matches("(dragmap)")) {
    BWADir = Channel.empty()
    
    if(params.dragmap_index_directory!=null) {
      DragmapDir = [[id:"dragmapindex"],file(params.dragmap_index_directory)];
      }
    else
      {
      DRAGMAP_INDEX(
        PREPARE_ONE_REFERENCE.out.fasta,
        [[id:"nomask"],[]]
        )
      versions = versions.mix(DRAGMAP_INDEX.out.versions)
      DragmapDir = DRAGMAP_INDEX.out.dragmap_index
      }
    }
   else
    {
    throw new IllegalArgumentException("unknown mapper: ${params.mapper}");
    }

  /***************************************************
   *
   *  MAP FASTQS TO BAM
   */
  if(params.mapper.matches("(pb|parabricks)")) {

    single_end
      .map{throw new IllegalArgumentException("SORRY: parabricks does'nt support single-end (2025-11-12)");}

    /***************************************************
    *
    *  MAP FASTQS TO BAM USING PARABRICKS
    *
    */
    PB_FQ2BAM(
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      BWADir,
      known_indels_ch,
      prevent_multi_gpu_ch //prevent multiple GPUs jobs at the same time
        .combine(paired_end)
        .map{it[1..<it.size()]} // remove item[0] containing multiple GPLUs guard
      )  
    versions = versions.mix(PB_FQ2BAM.out.versions)
    multiqc = multiqc.mix(PB_FQ2BAM.out.duplicate_metrics)
    bams_ch =  PB_FQ2BAM.out.bam
    prevent_multi_gpu_ch = bams_ch.count()
    }
  else if(params.mapper.matches("(bwa)")) {
    /***************************************************
    *
    *  MAP FASTQS TO BAM USING BWA
    *
    */
    MAP_BWA(
      metadata.plus(
        with_fastp : false, /* already done before */
        with_cram :  parseBoolean(params.with_cram),
        fastqc_before : parseBoolean(params.with_fastqc),
        fastqc_after : parseBoolean(params.with_fastqc),
        with_seqkit_split : parseBoolean(params.with_seqkit_split2),
        with_bqsr : parseBoolean(params.with_bqsr),
        markdup_enabled : parseBoolean(params.with_markdup),
        markdup_method : params.markdup_method
        ) ,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      BWADir,
      [[id:"bqsr_files"],[]],
      [[id:"nobed"],[]],
      single_end.mix(paired_end)
      )
    versions = versions.mix(MAP_BWA.out.versions)
    bams_ch =  MAP_BWA.out.bams
    }
  else if(params.mapper.matches("(dragmap)")) {
    MAP_DRAGMAP(
      metadata.plus(
            with_fastp : parseBoolean(params.with_fastp),
            with_cram :  parseBoolean(params.with_cram),
            fastqc_before : parseBoolean(params.with_fastqc),
            fastqc_after : parseBoolean(params.with_fastqc),
            with_seqkit_split : parseBoolean(params.with_seqkit_split2),
            with_bqsr : parseBoolean(params.with_bqsr),
            markdup_enabled : parseBoolean(params.with_markdup),
            markdup_method : params.markdup_method
            ) ,
          PREPARE_ONE_REFERENCE.out.fasta,
          PREPARE_ONE_REFERENCE.out.fai,
          PREPARE_ONE_REFERENCE.out.dict,
          DragmapDir,
          [[id:"bqsr_files"],[]],
          [[id:"nobed"],[]],
          single_end.mix(paired_end)
          )
    versions = versions.mix(MAP_DRAGMAP.out.versions)
    bams_ch =  MAP_DRAGMAP.out.bams
    }
  else
    {
      throw new IllegalArgumentException("unknown mapper: ${params.mapper}");
    }


  /* add extra BAM , if any */
  if(params.exta_bams!=null) {
      bams2_ch =  Channel.fromPath(params.exta_bams)
      if(params.exta_bams.endsWith(".json")) {
        bams2_ch = bams2_ch.splitJSon()
      } else {
        bams2_ch = bams2_ch.splitCsv(sep:',',header:true)
      }
      META_TO_BAMS(
        metadata ,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        bams2_ch
        )
      versions = versions.mix(META_TO_BAMS.out.versions)
      bams_ch = bams_ch.mix(META_TO_BAMS.out.bams)
      }

  /***************************************************
   *
   * Check no duplicate sample after merging
   *
   */  
  bams_ch.map{it[0].id}.unique().count()
    .combine( bams_ch.map{it[0].id}.count() )
    .map{
      if(!it[0].equals(it[1])) throw new IllegalArgumentException("probably duplicate sample. Check samplesheet");
      return it;
    }

  


  /***************************************************
   *
   * PREPARE THE REFERENCE GENOME
   *
   */
  if(params.bed==null) {
    if(is_exome) throw new IllegalArgumentException("EXOME but no capture was defined");
		bed = Channel.of([ [id:"nobed"], [] ])
		}
	else {
		PREPARE_USER_BED(
			metadata.plus([:]),
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			PREPARE_ONE_REFERENCE.out.scatter_bed,
			Channel.of([[id:"capture"],file(params.bed)])
			)
		versions = versions.mix(PREPARE_USER_BED.out.versions)
		multiqc = multiqc.mix(PREPARE_USER_BED.out.multiqc)
		bed = PREPARE_USER_BED.out.bed
		}
  
  /***************************************************
   *
   * STATS ON BAM
   *
   */
  // stats from mosdepth will be used by graphtyper
  mosdepth_summary_ch = Channel.empty()

  def run_bam_qc = parseBoolean(params.with_bam_qc) ||
        (params.mosdepth_min_depth as double) > 0.0 ||
        parseBoolean(params.with_graphtyper);

  // stats are required if graphtyper is used
  if( run_bam_qc ) {
      BAM_QC(
        metadata.plus([with_collect_metrics: parseBoolean(params.with_collect_metrics)]),
        PREPARE_ONE_REFERENCE.out.fasta,
			  PREPARE_ONE_REFERENCE.out.fai,
			  PREPARE_ONE_REFERENCE.out.dict,
        (params.bed==null?PREPARE_ONE_REFERENCE.out.scatter_bed:bed),
        bams_ch
        )
      versions = versions.mix(BAM_QC.out.versions)
      multiqc = multiqc.mix(BAM_QC.out.multiqc)
      mosdepth_summary_ch = BAM_QC.out.mosdepth_summary
  }
  /***************************************************
   *
   *  filter samples on quality ? mosdepth row 'total' must be greater than treshold.
   *
   */
  if( (params.mosdepth_min_depth as double) > 0.0 && run_bam_qc) {
    qual_bam_ch = mosdepth_summary_ch
        .splitCsv(header:true,sep:'\t')
        .filter{meta,row->row.chrom=="total"}
        .map{meta,row->meta.plus(coverage_depth: (row.chrom.mean as double)) }
        .map{meta->[meta, meta.coverage_depth >= (params.mosdepth_min_depth as double)] }
        .map{meta,depth_status->[meta.id,meta,depth_status]}
        .join(bams_ch.map{meta,bam,bai->[meta.id,meta,bam,bai]}, failOnMismatch:true, failOnDuplicate:true)
        .map{_sample_id,meta1,depth_status,meta2,_bam,bai->[meta1,depth_status,meta2,bam,bai]}
  
    // emit a warning if depth is too low
    qual_bam_ch
      .filter{_meta1,depth_status,_meta2,_bam,_bai->depth_status==false}
      .map{meta1,depth_status,meta2,bam,bai->log.warn("The following BAM have a low depth: ${meta1} ${bam}"); return 0;}

    if(params.on_low_bam_quality=="skip") {
      bams_ch = qual_bam_ch
        .filter{_meta1,depth_status,_meta2,_bam,_bai->depth_status==true}
        .map{_meta1,_depth_status,meta2,bam,bai->[meta2,bam,bai]}
      }
    else if(params.on_low_bam_quality=="abort") {
      bams_ch =  bams_ch .count()
            .combine(qual_bam_ch
                .filter{meta1,depth_status,meta2,bam,bai->depth_status==true}
                .count()
                )
            .filter{N1,N2->N1==N2} //same number of bam after filtering on DEPTH
            .combine(bams_ch)
            .map{_N1,_N2,meta,bam,bai->[meta,bam,bai]}
      }
    else
      {
      // do nothing
      }
  }


  /***************************************************
   *
   * UPDATE SEX to meta if unknown
   *
   */ 
  if( run_bam_qc) {
   bams_ch = bams_ch.map{meta,bam,bai->[meta.id,meta,bam,bai]}
      .join(BAM_QC.out.sample_sex_ch.map{meta->[meta.id,meta]}, failOnMismatch:true,failOnDuplicate:true )
      .map{_sample_id,meta1,bam,bai,meta2->
        if(!isBlank(meta1.sex) && meta1.sex.matches("(male|female)")) {
           if(!isBlank(meta2.sex) && meta2.sex.matches("(male|female)") && meta1.sex!=meta2.sex) {
            log.warn("inferred sex is not the same as declared sex ${meta1} ${meta2}");
            }
          return [meta1,bam,bai];
          }
        if(!isBlank(meta2.sex) && meta2.sex.matches("(male|female)")) return [meta1.plus(sex:meta2.sex),bam,bai];
        return [meta1,bam,bai];
      }
  }
  

  /***************************************************
   *
   * BUILD A PEDIGREE FROM META INFO
   *
   */
  META_TO_PED(metadata, bams_ch.map{it[0]})
	versions = versions.mix(META_TO_PED.out.versions)

  pedigree_ch = META_TO_PED.out.pedigree_gatk
  trios_ped_ch = META_TO_PED.out.pedigree_gatk //TODO check this, should be a pedigree of trios

  /***************************************************
   *
   * SOMALIER
   *
   */
  if(parseBoolean(params.with_somalier) && is_wgs) {
    // vcf for somalier
    if(params.somalier_vcf_sites!=null) {
        SOMALIER_SITES_VCF_INPUT(
          metadata.plus([
            arg_name : "somalier_vcf_sites",
            path: params.somalier_vcf_sites,
            require_index : true,
            required : true,
            unique : true
            ])
          )
        versions = versions.mix(SOMALIER_SITES_VCF_INPUT.out.versions)
        user_sites = SOMALIER_SITES_VCF_INPUT.out.vcf.first()
        }
      else
        {
        user_sites = Channel.of([[id:"no_user_sites_somalier"],[],[]]).first()
        }


    SOMALIER_BAMS(
          metadata,
          PREPARE_ONE_REFERENCE.out.fasta,
          PREPARE_ONE_REFERENCE.out.fai,
          PREPARE_ONE_REFERENCE.out.dict,
          bams_ch,
          META_TO_PED.out.pedigree_gatk,
          user_sites
          )
    versions = versions.mix(SOMALIER_BAMS.out.versions)
    multiqc = multiqc.mix(SOMALIER_BAMS.out.multiqc)
    }
  
  /***************************************************
   *
   * MANTA
   *
   */
  if(parseBoolean(params.with_manta)) {
    if(is_wgs) {
        manta_bed = [[id:"no_bed"],[]]
        }
    else
        {
        manta_bed =  bed
        }

    MANTA(
        metadata.plus(
          with_truvari : parseBoolean(params.with_truvari)
          ),
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
			  PREPARE_ONE_REFERENCE.out.dict,
        bams_ch,
        manta_bed
        )
    versions = versions.mix(MANTA.out.versions)
    multiqc = multiqc.mix(MANTA.out.multiqc)
    }

/***************************************************
   *
   * SMOOVE
   *
   */
  if(parseBoolean(params.with_smoove) && is_wgs) {
    SMOOVE(
        metadata,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        PREPARE_ONE_REFERENCE.out.complement_bed,
        bams_ch
        )
    versions = versions.mix(SMOOVE.out.versions)
    }

  /***************************************************
   *
   * INDEXCOV
   *
   */
  if(parseBoolean(params.with_indexcov) && is_wgs) {
    INDEXCOV(
      metadata.plus(batch_size: params.indexcov_batchsize),
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      PREPARE_ONE_REFERENCE.out.complement_bed,
      META_TO_PED.out.pedigree_gatk,
      bams_ch.map{meta,bam,bai->[meta,bam]/* no BAI please */}
      )
    versions = versions.mix(INDEXCOV.out.versions)
    multiqc = multiqc.mix(INDEXCOV.out.multiqc)
  }
  /***************************************************
   *
   * DELLY2
   *
   */
  if(parseBoolean(params.with_delly) && is_wgs) {
      DELLY(
            metadata,
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            PREPARE_ONE_REFERENCE.out.complement_bed,//TODO blacklisted region, use something better
            bams_ch
            )
      versions = versions.mix(DELLY.out.versions)
      multiqc = multiqc.mix(DELLY.out.multiqc)
  }


 /***************************************************
   *
   * CNVNATOR
   *
   */
  if(parseBoolean(params.with_cnvnator) && is_wgs) {
    CNVNATOR(
      metadata,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      bams_ch
      )
    versions = versions.mix(CNVNATOR.out.versions)
     multiqc = multiqc.mix(CNVNATOR.out.multiqc)
    }
 /***************************************************
   *
   * Download encode black list
   *
   */
  ENCODE_BLACKLIST(PREPARE_ONE_REFERENCE.out.dict)
  versions = versions.mix(ENCODE_BLACKLIST.out.versions)

  snv_calling_bed =  (params.bed!=null?bed:PREPARE_ONE_REFERENCE.out.scatter_bed)

 /***************************************************
   *
   * Remove blacklisted region for calling
   *
   */
  if(parseBoolean(params.exclude_encode_blacklist)) {
    BEDTOOLS_SUBTRACT(snv_calling_bed.combine(ENCODE_BLACKLIST.out.bed).map{meta1,bed1,meta2,bed2->[meta1,bed1,bed2]})
    versions = versions.mix(BEDTOOLS_SUBTRACT.out.versions)
    snv_calling_bed = BEDTOOLS_SUBTRACT.out.bed.first()
    }


  /* cut the bed/genome into parts for SNV calling per region */
  BEDTOOLS_MAKEWINDOWS( snv_calling_bed )
  versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

  /* if it's an exome , group the small genome together in BED */
  BED_CLUSTER1(
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      BEDTOOLS_MAKEWINDOWS.out.bed
      )
  versions = versions.mix(BED_CLUSTER1.out.versions)
  cluster_bed1 = BED_CLUSTER1.out.bed
    .map{meta,bed->bed}
    .map{it instanceof List?it:[it]}
    .flatMap()
    .map{bed->[[id:bed.baseName],bed]}

  /* But for slow callers like bcftools or freebayes, we need smaller intervals */
  BED_CLUSTER2(
    PREPARE_ONE_REFERENCE.out.fasta,
    PREPARE_ONE_REFERENCE.out.fai,
     PREPARE_ONE_REFERENCE.out.dict,
    BEDTOOLS_MAKEWINDOWS.out.bed
    )
  versions = versions.mix(BED_CLUSTER2.out.versions)
  cluster_bed2 = BED_CLUSTER2.out.bed.map{it[1]}
    .flatMap()
    .map{bed->[[id:bed.baseName],bed]}// give a uniq name to the beds
  
  grouped_bams = bams_ch
    .map{[it[1],it[2]]}//bam,bai
    .collect()
    .map{[[id:"bams"],it.flatten().sort()]}

  vcf_ch = Channel.empty()

 /***************************************************
   *
   * PARABRICKS HAPLOTYPECALLER
   *
   */
  if(parseBoolean(params.with_pb_haplotypecaller)) {
    PB_HAPLOTYPECALLER(
      metadata,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      cluster_bed1, // not for calling but joining
      [[id:"no_glnexus_config"],[]],
      prevent_multi_gpu_ch //prevent multiple GPUs jobs at the same time
        .combine(  bams_ch.map{it + [[]]/* empty file for BQSR */})
        .map{it[1..<it.size()]} // remove item[0] containing multiple GPLUs guard
      )
    versions = versions.mix(PB_HAPLOTYPECALLER.out.versions)
    vcf_ch = vcf_ch.mix(PB_HAPLOTYPECALLER.out.vcf)
    prevent_multi_gpu_ch = vcf_ch.count()
  }

 /***************************************************
   *
   * PARABRICKS DEEPVARIANT
   *
   */
  if(parseBoolean(params.with_pb_deepvariant)) {
    PB_DEEPVARIANT(
      metadata,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      cluster_bed1, // not for calling but joining
      prevent_multi_gpu_ch //prevent multiple GPUs jobs at the same time
        .combine(bams_ch)
        .map{it[1..<it.size()]} // remove item[0] containing multiple GPLUs guard
      )
    versions = versions.mix(PB_DEEPVARIANT.out.versions)
    vcf_ch = vcf_ch.mix(PB_DEEPVARIANT.out.vcf)
    prevent_multi_gpu_ch =  vcf_ch.count()
  }


  /***************************************************
   *
   *  GRAPHTYPER
   *
   */
  if(parseBoolean(params.with_graphtyper)) {
    GRAPHTYPER(
      metadata,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      mosdepth_summary_ch,
      cluster_bed2,
      bams_ch
      )
    versions = versions.mix(GRAPHTYPER.out.versions)
    vcf_ch = vcf_ch.mix(GRAPHTYPER.out.vcf)
  }

  /***************************************************
   *
   *  BCFTOOLS mpileup/call
   *
   */
  if(parseBoolean(params.with_bcftools_call)) {
    BCFTOOLS_CALL(
      metadata,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      pedigree_ch,
      cluster_bed2,
      bams_ch
      )
    versions = versions.mix(BCFTOOLS_CALL.out.versions)
    vcf_ch = vcf_ch.mix(BCFTOOLS_CALL.out.vcf)
  }

/***************************************************
   *
   *  FREEBAYES call
   *
   */
  if(parseBoolean(params.with_freebayes)) {
    FREEBAYES_CALL(
      metadata,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      pedigree_ch,
      cluster_bed2,
      bams_ch
      )
    versions = versions.mix(FREEBAYES_CALL.out.versions)
    vcf_ch = vcf_ch.mix(FREEBAYES_CALL.out.vcf)
    }

  if(parseBoolean(params.with_haplotype_caller)) {
    HAPLOTYPECALLER(
      metadata.plus(gvcf_merge_method:params.hc_gvcf_merge_method),
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      [[id:"refs"],[]],//other ref TODO
      dbsnp_ch,
      [[id:"ped"],[]], //TODO ped
      cluster_bed1,
      bams_ch
      )
    versions = versions.mix(HAPLOTYPECALLER.out.versions)
    vcf_ch = vcf_ch.mix(HAPLOTYPECALLER.out.vcf)
    }   

  /**
   *
   * UNMAPPED reads with spades
   *
   */
  if(parseBoolean(params.with_unmapped)) {
    UNMAPPED(
      metadata,
      [[id:"contaminants"],file("${moduleDir}/../unmapped/contaminants.txt")],
      bams_ch
        .combine(PREPARE_ONE_REFERENCE.out.fasta)
        .combine(PREPARE_ONE_REFERENCE.out.fai)
        .combine(PREPARE_ONE_REFERENCE.out.dict)
        .map{meta1,bam,bai,m2,fa,m3,fai,m4,dict->[meta1,bam,bai,fa,fai,dict]}
      )
    versions = versions.mix(UNMAPPED.out.versions)
    multiqc = multiqc.mix(UNMAPPED.out.multiqc)
  }

  /***************************************************
   *
   * SKIP EMPTY VCFS for stats
   *
   */
  NOT_EMPTY_VCF(vcf_ch)
  versions = versions.mix(NOT_EMPTY_VCF.out.versions)
  vcf_ch = NOT_EMPTY_VCF.out.vcf
    .filter{meta,vcf,idx,status->status.contains("VARIANT")}
    .map{meta,vcf,idx,status->[meta,vcf,idx]}


  /***************************************************
   *
   * GUESS PLOIDY FROM VCFS
   *
   */
  BCFTOOLS_GUESS_PLOIDY(
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      vcf_ch
      )
  versions = versions.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)

  /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
  BCFTOOLS_STATS(
    PREPARE_ONE_REFERENCE.out.fasta,
    PREPARE_ONE_REFERENCE.out.fai,
    bed,
    gtf_ch.map{meta,gtf,idx->[meta,gtf]}.first(),//meta,gtf
    [[id:"no_samples"],[]],
    vcf_ch.map{[it[0],[it[1],it[2]]]}
    )
  versions = versions.mix(BCFTOOLS_STATS.out.versions)

  /***************************************************
   *
   *  MULTIQC
   *
   */
  if( parseBoolean(params.with_multiqc)) {
    COMPILE_VERSIONS(versions.collect().map{it.sort()})
    multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MULTIQC([[id:"no_mqc_config"],[]],
      multiqc.map{it[1]}.collect().map{[[id:"parabricks"],it]})
    }
}

runOnComplete(workflow)

/***************************************************
  *
  *  BUILD A PEDIGREE FILE
  *
  */
process MAKE_PED_xxxxxxxxxxxx {
executor "local"
input:
  tuple val(meta),val(content)
output:
  tuple val(meta),path("pedigree.ped"),emit:pedigree
  tuple val(meta),path("trios.ped"),optional:true,emit:trio_ped
  tuple val(meta),path("gatk.ped"),optional:true,emit:gatk_ped

script:
"""
mkdir -p TMP

cat << __EOF__  > pedigree.ped
${content}
__EOF__

# special pedigree containing at least one trios
awk '(!(\$3=="0" || \$4=="0"))' pedigree.ped > TMP/a

if test -s TMP/jeter.a
then
  cp TMP/jeter.ped trios.ped
fi

awk -f "${moduleDir}/../../modules/gatk/possibledenovo/pedigree4gatk.awk" pedigree.ped > gatk.ped

"""
}
