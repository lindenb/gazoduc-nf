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
include {PB_HAPLOTYPECALLER                       } from '../../modules/parabricks/haplotypecaller'
include {SCATTER_TO_BED                           } from '../../subworkflows/gatk/scatterintervals2bed'
include {BAM_QC                                   } from '../../subworkflows/bamqc'
include {SOMALIER_BAMS                            } from '../../subworkflows/somalier/bams'
include {MANTA_SINGLE                             } from '../../subworkflows/manta/single'
include {INDEXCOV                                 } from '../../subworkflows/indexcov/simple'
include {DELLY                                    } from '../../subworkflows/delly2'
include {CNVNATOR                                 } from '../../subworkflows/cnvnator'
include {BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include {BED_CLUSTER                              } from '../../modules/jvarkit/bedcluster'
include {GRAPHTYPER                               } from '../../modules/graphtyper/genotype'
include {BCFTOOLS_CALL                            } from '../../subworkflows/bcftools/call'
include {FREEBAYES_CALL                           } from '../../subworkflows/freebayes/call'


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
      BAM_QC(hash_ref,fasta,fai,dict,bed,PB_FQ2BAM.out.bam)
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
            PB_FQ2BAM.out.bam,
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
        PB_FQ2BAM.out.bam,
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
            PB_FQ2BAM.out.bam.map{[it[0],it[1]]/* no BAI please */}
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
            PB_FQ2BAM.out.bam
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
      PB_FQ2BAM.out.bam
      )
    versions_ch = versions_ch.mix(CNVNATOR.out.versions)
    }

 /***************************************************
   *
   * HAPLOTYPECALLER
   *
   */
  if(params.with_haplotypecaller == true) {
    PB_HAPLOTYPECALLER(
      fasta,
      fai,
      PB_FQ2BAM.out.bam.map{it + [[]]/* empty file for BQSR */}
      )
    versions_ch = versions_ch.mix(PB_HAPLOTYPECALLER.out.versions)
  }

  /* cut the bed/genome into parts for SV calling per region */
  BEDTOOLS_MAKEWINDOWS(bed)
  versions_ch = versions_ch.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

  /* if it's an exome , group the small genome together in BED */
  BED_CLUSTER(fasta,fai,dict,BEDTOOLS_MAKEWINDOWS.out.bed)
  versions_ch = versions_ch.mix(BED_CLUSTER.out.versions)
  cluster_bed = BED_CLUSTER.out.bed.map{it[1]}.flatMap()

  grouped_bams = PB_FQ2BAM.out.bam
    .map{[it[1],it[2]]}//bam,bai
    .collect()
    .map{[[id:"bams"],it.flatten().sort()]}

  /***************************************************
   *
   *  GRAPHTYPER
   *
   */
  if(params.with_graphtyper == true) {
    GRAPHTYPER(
      fasta,
      fai,
      grouped_bams.combine(cluster_bed).view()
      )
    versions_ch = versions_ch.mix(GRAPHTYPER.out.versions)
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
      PB_FQ2BAM.out.bam
      )
    versions_ch = versions_ch.mix(BCFTOOLS_CALL.out.versions)
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
      PB_FQ2BAM.out.bam
      )
    versions_ch = versions_ch.mix(FREEBAYES_CALL.out.versions)
    }

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
