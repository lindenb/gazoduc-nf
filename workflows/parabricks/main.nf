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
include {FASTQC as FASTQC_BEFORE_TRIM             } from '../../modules/fastqc'
include {FASTQC as FASTQC_AFTER_TRIM              } from '../../modules/fastqc'
include {FASTP                                    } from '../../modules/fastp'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions'
include {BWA_INDEX                                } from '../../modules/bwa/index'
include {PB_FQ2BAM                                } from '../../modules/parabricks/fq2bam'
include {PB_HAPLOTYPECALLER                       } from '../../subworkflows/parabricks/haplotypecaller'
include {PB_DEEPVARIANT                           } from '../../subworkflows/parabricks/deepvariant'
include {SCATTER_TO_BED                           } from '../../subworkflows/gatk/scatterintervals2bed'
include {BAM_QC                                   } from '../../subworkflows/bamqc'
include {SOMALIER_BAMS                            } from '../../subworkflows/somalier/bams'
include {MANTA_SINGLE                             } from '../../subworkflows/manta/single'
include {INDEXCOV                                 } from '../../subworkflows/indexcov/simple'
include {DELLY                                    } from '../../subworkflows/delly2'
include {CNVNATOR                                 } from '../../subworkflows/cnvnator'
include {BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include {BED_CLUSTER                              } from '../../modules/jvarkit/bedcluster'
include {BCFTOOLS_GUESS_PLOIDY                    } from '../../modules/bcftools/guess_ploidy'
include {BCFTOOLS_STATS                           } from '../../modules/bcftools/stats'
include {GRAPHTYPER                               } from '../../modules/graphtyper/genotype'
include {BCFTOOLS_CALL                            } from '../../subworkflows/bcftools/call'
include {FREEBAYES_CALL                           } from '../../subworkflows/freebayes/call'
include {HAPLOTYPECALLER                          } from '../../subworkflows/gatk/haplotypecaller'

Map assertKeyExists(final Map hash,final String key) {
    if(!hash.containsKey(key)) throw new IllegalArgumentException("no key ${key}'in ${hash}");
    return hash;
}

Map assertKeyExistsAndNotEmpty(final Map hash,final String key) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(value.isEmpty()) throw new IllegalArgumentException("empty ${key}'in ${hash}");
    return hash;
}

Map assertKeyMatchRegex(final Map hash,final String key,final String regex) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(!value.matches(regex)) throw new IllegalArgumentException(" ${key}'in ${hash} doesn't match regex '${regex}'.");
    return hash;
}


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {

  def hash_ref= [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
                  ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ hash_ref, file(params.fasta)]
	def fai   = [ hash_ref, file(params.fai)]
	def dict  = [ hash_ref, file(params.dict)]
  def known_indels = [hash_ref, [],[]]
  def gtf     = [hash_ref,[],[]]
  
  if(params.gtf!=null) {
    gtf = [hash_ref,file(params.gtf),file(params.gtf+".tbi")]
  }


  if(!(params.hts_type.equals("WGS") || params.hts_type.equals("WES"))) {
    throw new IllegalArgumentException("hts_type should batch WES|WGS ${params.hts_type}")
  }
  def is_exome = params.hts_type.equals("WES");
  def is_wgs = params.hts_type.equals("WGS");

  versions_ch = Channel.empty()
  multiqc_ch = Channel.empty()

  	samplesheet_ch = Channel.fromPath(params.samplesheet)
			.splitCsv(sep:',',header:true)
			.map{assertKeyMatchRegex(it,"sample","^[A-Za-z][A-Za-z0-9_\\-\\.]*[A-Za-z0-9]\$")}
			.map{assertKeyMatchRegex(it,"fastq_1","^\\S+\\.(fastq|fq|ora|fastq\\.gz|fq\\.gz)\$")}
      .map{
          if(it.containsKey("fastq_2")) {
            if(!(it.fastq_2==null || it.fastq_2.equals(".") || it.fastq_2.isEmpty())) {
              // fastq cannot end with '.ora'
              assertKeyMatchRegex(it,"fastq_2","^\\S+\\.(fastq|fq|fastq\\.gz|fq\\.gz)\$")
              }
            return it;
            }
          }
		.map{[ [id:it.sample],it.fastq_1,it.fastq_2]}
			

	if(params.sampleinfo!=null) {
		ch2 = Channel.fromPath(params.sampleinfo)
                        .splitCsv(sep:',',header:true)
			.map{assertKeyExistsAndNotEmpty(it,"sample")}
			.map{[ [id:it.sample], it ]}

		/* join with fastq info */
		samplesheet_ch = samplesheet_ch.join(ch2,remainder:true)
			.filter{it.size()==4}
			.map {
			  if(it[3]!=null) return [it[0].plus(it[3]),it[1],it[2]];
			  return  [it[0],it[1],it[2]]	
         }
		}
	   


    /* select fastq having an .ORA extension */
    fastq_type_ch = samplesheet_ch.branch{v->
      ora       : v[1].endsWith(".ora")
      single_end: v[2]==null || v[2].isEmpty() || v[2].equals(".")
      paired_end: true
      }

  /***************************************************
   *
   *  CONVERT ORA TO FASTQ
   *
   */
  ORA_TO_FASTQ(
    Channel.of([id:"ora"]),
    fastq_type_ch.ora.map{[it[0],file(it[1])]}
    ) 
  versions_ch = versions_ch.mix(ORA_TO_FASTQ.out.versions)    

  single_end = fastq_type_ch.single_end
		.map{[it[0], file(it[1])]}
		.mix(ORA_TO_FASTQ.out.single_end)
    .map{[it[0],[it[1]]]}
    
    
  paired_end = fastq_type_ch.paired_end
		.mix(ORA_TO_FASTQ.out.paired_end)
    .map{[it[0],[it[1],it[2]]]}

   
  /***************************************************
   *
   *  FASTQC : get FASTQ quality control
   *
   */
   if( params.with_fastqc == true) {
      FASTQC_BEFORE_TRIM( single_end.mix(paired_end)
      )
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

  /***************************************************
   *
   *  Index FASTA For BWA
   *
   */
  BWA_INDEX(fasta)
  versions_ch = versions_ch.mix(BWA_INDEX.out.versions)
  
  /***************************************************
   *
   *  MAP FASTQS TO BAM USING PARABRICKS
   *
   */
  PB_FQ2BAM(
    fasta,
    fai,
    BWA_INDEX.out.bwa_index ,
    known_indels,
    single_end.mix(paired_end)
  )  
  versions_ch = versions_ch.mix(PB_FQ2BAM.out.versions)
  multiqc_ch = multiqc_ch.mix(PB_FQ2BAM.out.duplicate_metrics)

  bams_ch =  PB_FQ2BAM.out.bam


  /* add extra BAM , if any */
  if(params.exta_bams!=null) {
      bams2_ch = 
        Channel.fromPath(params.exta_bams)
          .splitCsv(sep:',',header:true)
          .map{assertKeyMatchRegex(it,"sample","^[A-Za-z][A-Za-z0-9_\\-\\.]*[A-Za-z0-9]\$")}
          .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
          .map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
          .map{[ [id:it.sample],it.bam,it.bai]}
      
        /* join SAMPLE INFO is any */
      	if(params.sampleinfo!=null) {
            ch2 = Channel.fromPath(params.sampleinfo)
                                .splitCsv(sep:',',header:true)
              .map{assertKeyExistsAndNotEmpty(it,"sample")}
              .map{[ [id:it.sample], it ]}

            /* join with fastq info */
            bams2_ch = bams2_ch.join(ch2,remainder:true)
              .filter{it.size()==4}
              .map {
                if(it[3]!=null) return [it[0].plus(it[3]),it[1],it[2]];
                return  [it[0],it[1],it[2]]	
                }
            }


       bams_ch = bams_ch.mix(bams2_ch)
  }

  /***************************************************
   *
   * GET A BED FOR GENOME WITHOUT THE POLY-NNNN
   * will be used for bam QC
   *
   */
  if(params.bed==null) {
    SCATTER_TO_BED(hash_ref,fasta,fai,dict)
    versions_ch = versions_ch.mix(SCATTER_TO_BED.out.versions)
    bed = SCATTER_TO_BED.out.bed
  } else {
      bed = [hash_ref, file(params.bed)]
  }
  /***************************************************
   *
   * STATS ON BAM
   *
   */
  if( params.with_bam_qc == true) {
      BAM_QC(hash_ref,fasta,fai,dict,bed,bams_ch)
      versions_ch = versions_ch.mix(BAM_QC.out.versions)
      multiqc_ch = multiqc_ch.mix(BAM_QC.out.multiqc)
  }




  /***************************************************
   *
   * SOMALIER
   *
   */ 
  if(params.with_somalier==true && is_wgs) {
        pedigree = [[:],[]] //TODO check pedigree
        user_sites = [[:],[],[]] //custom sites
        SOMALIER_BAMS(
            hash_ref,
            Channel.of(fasta),
            fai,
            dict,
            bams_ch,
            pedigree,
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
        fasta,
        fai,
        dict,
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
     pedigree = [[:],[]] //TODO check pedigree
    INDEXCOV(
            hash_ref,
            fasta,
            fai,
            dict,
            pedigree,
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
            fasta,
            fai,
            dict,
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
      fasta,
      fai,
      dict,
      bams_ch
      )
    versions_ch = versions_ch.mix(CNVNATOR.out.versions)
    }



  /* cut the bed/genome into parts for SV calling per region */
  BEDTOOLS_MAKEWINDOWS(bed)
  versions_ch = versions_ch.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

  /* if it's an exome , group the small genome together in BED */
  BED_CLUSTER(fasta,fai,dict,BEDTOOLS_MAKEWINDOWS.out.bed)
  versions_ch = versions_ch.mix(BED_CLUSTER.out.versions)
  cluster_bed = BED_CLUSTER.out.bed.map{it[1]}.flatMap()

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
      fasta,
      fai,
      dict,
      cluster_bed.map{[[id:it.name.toString().md5().substring(0,8)],it]},
      bams_ch.map{it + [[]]/* empty file for BQSR */}
      )
    versions_ch = versions_ch.mix(PB_HAPLOTYPECALLER.out.versions)
    vcf_ch = vcf_ch.mix(PB_HAPLOTYPECALLER.out.vcf)
  }

 /***************************************************
   *
   * PARABRICKS DEEPVARIANT
   *
   */
  if(params.with_pb_deepvariant == true) {
    PB_DEEPVARIANT(
      hash_ref,
      fasta,
      fai,
      dict,
      cluster_bed.map{[[id:it.name.toString().md5().substring(0,8)],it]},
      bams_ch
      )
    versions_ch = versions_ch.mix(PB_DEEPVARIANT.out.versions)
    vcf_ch = vcf_ch.mix(PB_DEEPVARIANT.out.vcf)
  }


  /***************************************************
   *
   *  GRAPHTYPER
   *
   */
  if(params.with_graphtyper == true) {
    GRAPHTYPER(
      fasta,
      fai,
      grouped_bams.combine(cluster_bed)
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
    pedigree = [[:],[]] //TODO check pedigree
    BCFTOOLS_CALL(
      hash_ref,
      fasta,
      fai,
      dict,
      pedigree,
      cluster_bed.map{[[id:"bed"],it]},
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
    pedigree = [[:],[]] //TODO check pedigree
    FREEBAYES_CALL(
      hash_ref,
      fasta,
      fai,
      dict,
      pedigree,
      cluster_bed.map{[[id:"bed"],it]},
      bams_ch
      )
    versions_ch = versions_ch.mix(FREEBAYES_CALL.out.versions)
    vcf_ch = vcf_ch.mix(FREEBAYES_CALL.out.vcf)
    }

  if(params.with_gatk == true) {
    pedigree = [[:],[]] //TODO check pedigree
    HAPLOTYPECALLER(
      hash_ref,
      fasta,
      fai,
      dict,
      cluster_bed.map{[[id:"bed"],it]},
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
  BCFTOOLS_GUESS_PLOIDY(fasta, fai,vcf_ch)
  versions_ch = versions_ch.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)

  /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
  BCFTOOLS_STATS(
    fasta,
    fai,
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
    multiqc_ch = multiqc_ch.mix(COMPILE_VERSIONS.out.multiqc)
    multiqc_ch.view()
    //MULTIQC(multiqc_ch.collect().map{[[id:"parabricks"],it]})
    }
}

runOnComplete(workflow)
