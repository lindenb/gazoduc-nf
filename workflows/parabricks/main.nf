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
include {runOnComplete; dumpParams                } from '../../modules/utils/functions.nf'
include {ORA_TO_FASTQ                             } from '../../subworkflows/ora/ora2fastq'
include {FASTP                                    } from '../../subworkflows/fastp'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions'
include {BWA_INDEX                                } from '../../modules/bwa/index'
include {PB_FQ2BAM                                } from '../../modules/parabricks/fq2bam'
include {PB_HAPLOTYPECALLER                       } from '../../subworkflows/parabricks/haplotypecaller'
include {PB_DEEPVARIANT                           } from '../../subworkflows/parabricks/deepvariant'
include {PB_BAM2FQ                                } from '../../modules/parabricks/bam2fq'
include {SCATTER_TO_BED                           } from '../../subworkflows/gatk/scatterintervals2bed'
include {BAM_QC                                   } from '../../subworkflows/bamqc'
include {SOMALIER_BAMS                            } from '../../subworkflows/somalier/bams'
include {MANTA_SINGLE                             } from '../../subworkflows/manta/single'
include {INDEXCOV                                 } from '../../subworkflows/indexcov/simple'
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


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
  if(params.fasta==null) {
      throw new IllegalArgumentException("undefined --fasta");
      }
	
  if(!(params.hts_type.equals("WGS") || params.hts_type.equals("WES"))) {
    throw new IllegalArgumentException("hts_type should batch WES|WGS ${params.hts_type}")
    }
  
  def hash_ref= [
      id: file(params.fasta).baseName
      ]
  def metadata = hash_ref;
  def known_indels = [hash_ref, [],[]]
  def gtf     = [hash_ref,[],[]]
  
  if(params.gtf!=null) {
    gtf = [hash_ref,file(params.gtf),file(params.gtf+".tbi")]
  }


 
  def is_exome = params.hts_type.equals("WES");
  def is_wgs = params.hts_type.equals("WGS");

  versions_ch = Channel.empty()
  multiqc_ch = Channel.empty()

  /* no fastq samplesheet */
  if(params.samplesheet==null) {
    samplesheet0_ch = Channel.empty()
    }
  else
    {
    samplesheet0_ch = Channel.fromPath(params.samplesheet);
    if(params.samplesheet.endsWith(".json")) {
      samplesheet0_ch =  samplesheet0_ch.splitJson();
      }
    else
      {
      samplesheet0_ch =  samplesheet0_ch.splitCsv(sep:',',header:true);
      }
    }
  
  SAMPLESHEET_TO_FASTQ(
    metadata,
    samplesheet0_ch
    )
  versions_ch =  SAMPLESHEET_TO_FASTQ.out.versions


  single_end = SAMPLESHEET_TO_FASTQ.out.single_end
    .map{[it[0],[it[1]]]}
    
    
  paired_end = SAMPLESHEET_TO_FASTQ.out.paired_end
    .map{[it[0],[it[1],it[2]]]}

   
  /***************************************************
   *
   *  FASTQC : get FASTQ quality control
   *
   */
   if( params.with_fastqc == true) {
    FASTQC_BEFORE_TRIM( single_end.mix(paired_end) )
    versions_ch = versions_ch.mix(FASTQC_BEFORE_TRIM.out.versions)
    multiqc_ch = multiqc_ch.mix(FASTQC_BEFORE_TRIM.out.zip)
    }

  
  /***************************************************
   *
   *  FASTP trim treads
   *
   */
  if( params.with_fastp == true) {
    FASTP(
      single_end .mix(paired_end)
      )
    versions_ch = versions_ch.mix(FASTP.out.versions)
    multiqc_ch = multiqc_ch.mix(FASTP.out.json)


   /***************************************************
    *
    *  FASTQC : get FASTQ quality control after FASTP
    *
    */
    if( params.with_fastqc == true) {
      FASTQC_AFTER_TRIM( FASTP.out.fastqs)
      versions_ch = versions_ch.mix(FASTQC_AFTER_TRIM.out.versions)
      multiqc_ch = multiqc_ch.mix(FASTQC_AFTER_TRIM.out.zip)
      }
    // split single and paired ends 
    trim_reads = FASTP.out.fastqs.branch{v->
      paired_end : v[1].size()==2
      single_end: v[1].size()==1
      other: true
      }

    single_end = trim_reads.single_end
    paired_end = trim_reads.paired_end.map{[it[0],it[1].sort()]}
  }


	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta)
			.map{[[id:it.baseName],it]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


  /***************************************************
   *
   *  Index FASTA For BWA
   *
   */
  	if(params.bwa_index_directory!=null) {
      BWADir = [[id:"bwaindex"],file(params.bwa_index_directory)];
      }
    else
      {
      BWA_INDEX(PREPARE_ONE_REFERENCE.out.fasta)
      versions_ch = versions_ch.mix(BWA_INDEX.out.versions)
      BWADir = BWA_INDEX.out.bwa_index
      }
  /***************************************************
   *
   *  MAP FASTQS TO BAM
   */
  if(params.mapper.matches("(pb|parabricks)")) {
    /***************************************************
    *
    *  MAP FASTQS TO BAM USING PARABRICKS
    *
    */
    PB_FQ2BAM(
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      BWADir,
      known_indels,
      prevent_multi_gpu_ch //prevent multiple GPUs jobs at the same time
        .combine(single_end.mix(paired_end))
        .map{it[1..<it.size()]} // remove item[0] containing multiple GPLUs guard
      )  
    versions_ch = versions_ch.mix(PB_FQ2BAM.out.versions)
    multiqc_ch = multiqc_ch.mix(PB_FQ2BAM.out.duplicate_metrics)
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
      metadata ,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      BWADir,
      [[id:"bqsr_files"],[]],
      [[id:"nobed"],[]],
      single_end.mix(paired_end)
      )
  versions_ch = versions_ch.mix(MAP_BWA.out.versions)
  bams_ch =  MAP_BWA.out.bams
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
      versions_ch = versions_ch.mix(META_TO_BAMS.out.versions)
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
   * BUILD A PEDIGREE FROM META INFO
   *
   */
  META_TO_PED(metadata, bams_ch.map{it[0]})
	versions = versions.mix(META_TO_PED.out.versions)

   pedigree_ch = META_TO_PED.out.pedigree.ifEmpty([[:],[]])
   trios_ped_ch = META_TO_PED.out.trio_ped.ifEmpty([[:],[]])

  /***************************************************
   *
   * GET A BED FOR GENOME WITHOUT THE POLY-NNNN
   * will be used for bam QC
   *
   */
  if(params.bed==null) {
    if(is_exome) throw new IllegalArgimentException("EXOME but no capture was defined");
		bed = Channel.of([ [id:"nobed"], [] ])
		}
	else {
		PREPARE_USER_BED(
			metadata,
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
  // stats are required if graphtyper is used
  if( params.with_bam_qc == true || params.with_graphtyper == true) {
      BAM_QC(
        PREPARE_ONE_REFERENCE.out.fasta,
			  PREPARE_ONE_REFERENCE.out.fai,
			  PREPARE_ONE_REFERENCE.out.dict,
        bed,
        bams_ch
        )
      versions_ch = versions_ch.mix(BAM_QC.out.versions)
      multiqc_ch = multiqc_ch.mix(BAM_QC.out.multiqc)
      mosdepth_summary_ch = BAM_QC.out.mosdepth_summary
  }




  /***************************************************
   *
   * SOMALIER
   *
   */ 
  if(params.with_somalier==true && is_wgs) {
        user_sites = Channel.of([[:],[],[]]) //custom sites
        SOMALIER_BAMS(
            hash_ref,
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
			      PREPARE_ONE_REFERENCE.out.dict,
            bams_ch,
            META_TO_PED.out.pedigree,
            user_sites
            )
        versions_ch = versions_ch.mix(SOMALIER_BAMS.out.versions)
    }
  /***************************************************
   *
   * MANTA
   *
   */
  if(params.with_manta == true && is_wgs) {
    MANTA_SINGLE(
        hash_ref,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
			  PREPARE_ONE_REFERENCE.out.dict,
        bams_ch,
        [[:],[]] //No bed
        )
    versions_ch = versions_ch.mix(MANTA_SINGLE.out.versions)
    }

  /***************************************************
   *
   * INDEXCOV
   *
   */
  if(params.with_indexcov == true && is_wgs) {
    INDEXCOV(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      META_TO_PED.out.pedigree,
      bams_ch.map{[it[0],it[1]]/* no BAI please */}
      )
    versions_ch = versions_ch.mix(INDEXCOV.out.versions)
  }
  /***************************************************
   *
   * DELLY2
   *
   */
  if(params.with_delly == true && is_wgs) {
      DELLY(
            hash_ref,
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            bams_ch
            )
      versions_ch = versions_ch.mix(DELLY.out.versions)
  }


 /***************************************************
   *
   * CNVNATOR
   *
   */
  if(params.with_cnvnator == true) {
    CNVNATOR(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      bams_ch
      )
    versions_ch = versions_ch.mix(CNVNATOR.out.versions)
    }



  /* cut the bed/genome into parts for SV calling per region */
  BEDTOOLS_MAKEWINDOWS(bed)
  versions_ch = versions_ch.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

  /* if it's an exome , group the small genome together in BED */
  BED_CLUSTER1(fasta,fai,dict,BEDTOOLS_MAKEWINDOWS.out.bed)
  versions_ch = versions_ch.mix(BED_CLUSTER1.out.versions)
  cluster_bed1 = BED_CLUSTER1.out.bed.map{it[1]}.flatMap()

  /* But for slow callers like bcftools or freebayes, we need smaller intervals */
  BED_CLUSTER2(fasta,fai,dict,BEDTOOLS_MAKEWINDOWS.out.bed)
  versions_ch = versions_ch.mix(BED_CLUSTER2.out.versions)
  cluster_bed2 = BED_CLUSTER2.out.bed.map{it[1]}.flatMap()


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
  if(params.with_pb_haplotypecaller == true) {
    PB_HAPLOTYPECALLER(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      cluster_bed1.map{[[id:it.name.toString().md5().substring(0,8)],it]},
      prevent_multi_gpu_ch //prevent multiple GPUs jobs at the same time
        .combine(  bams_ch.map{it + [[]]/* empty file for BQSR */})
        .map{it[1..<it.size()]} // remove item[0] containing multiple GPLUs guard
      )
    versions_ch = versions_ch.mix(PB_HAPLOTYPECALLER.out.versions)
    vcf_ch = vcf_ch.mix(PB_HAPLOTYPECALLER.out.vcf)
    prevent_multi_gpu_ch = vcf_ch.count()
  }

 /***************************************************
   *
   * PARABRICKS DEEPVARIANT
   *
   */
  if(params.with_pb_deepvariant == true) {
    PB_DEEPVARIANT(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      cluster_bed1.map{[[id:it.name.toString().md5().substring(0,8)],it]},
      prevent_multi_gpu_ch //prevent multiple GPUs jobs at the same time
        .combine(bams_ch)
        .map{it[1..<it.size()]} // remove item[0] containing multiple GPLUs guard
      )
    versions_ch = versions_ch.mix(PB_DEEPVARIANT.out.versions)
    vcf_ch = vcf_ch.mix(PB_DEEPVARIANT.out.vcf)
    prevent_multi_gpu_ch =  vcf_ch.count()
  }


  /***************************************************
   *
   *  GRAPHTYPER
   *
   */
  if(params.with_graphtyper == true) {
    GRAPHTYPER(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      mosdepth_summary_ch,
      cluster_bed1.map{[[id:it.name.toString().md5().substring(0,8)],it]},
      bams_ch
      )
    versions_ch = versions_ch.mix(GRAPHTYPER.out.versions)
    vcf_ch = vcf_ch.mix(GRAPHTYPER.out.vcf)
  }

  /***************************************************
   *
   *  BCFTOOLS mpileup/call
   *
   */
  if(params.with_bcftools_call == true) {
    BCFTOOLS_CALL(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      pedigree_ch,
      cluster_bed2.map{[[id:"bed"],it]},
      bams_ch
      )
    versions_ch = versions_ch.mix(BCFTOOLS_CALL.out.versions)
    vcf_ch = vcf_ch.mix(BCFTOOLS_CALL.out.vcf)
  }

/***************************************************
   *
   *  FREEBAYES call
   *
   */
  if(params.with_freebayes == true) {
    FREEBAYES_CALL(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      pedigree_ch,
      cluster_bed2.map{[[id:"bed"],it]},
      bams_ch
      )
    versions_ch = versions_ch.mix(FREEBAYES_CALL.out.versions)
    vcf_ch = vcf_ch.mix(FREEBAYES_CALL.out.vcf)
    }

  if(params.with_gatk == true) {
    HAPLOTYPECALLER(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      PREPARE_ONE_REFERENCE.out.dict,
      cluster_bed1.map{[[id:"bed"],it]},
      bams_ch
      )
    versions_ch = versions_ch.mix(HAPLOTYPECALLER.out.versions)
    vcf_ch = vcf_ch.mix(HAPLOTYPECALLER.out.vcf)
    }   

  /***************************************************
   *
   * GUESS PLOIDY FROM VCFS
   *
   */
  BCFTOOLS_GUESS_PLOIDY(
      hash_ref,
      PREPARE_ONE_REFERENCE.out.fasta,
      PREPARE_ONE_REFERENCE.out.fai,
      vcf_ch
      )
  versions_ch = versions_ch.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)

  /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
  BCFTOOLS_STATS(
    PREPARE_ONE_REFERENCE.out.fasta,
    PREPARE_ONE_REFERENCE.out.fai,
    bed,
    Channel.of(gtf).map{[it[0],it[1]]}.first(),//meta,gtf
    [[:],[]],//samples,
    vcf_ch.map{[it[0],[it[1],it[2]]]}
    )
  versions_ch = versions_ch.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)

  /***************************************************
   *
   *  MULTIQC
   *
   */
  if( params.with_multiqc ==true) {
    COMPILE_VERSIONS(versions_ch.collect().map{it.sort()})
    multiqc_ch = multiqc_ch.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MULTIQC([[id:"no_mqc_config"],[]],
      multiqc_ch.map{it[1]}.collect().map{[[id:"parabricks"],it]})
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
