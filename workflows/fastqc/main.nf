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

include {dumpParams;runOnComplete   } from '../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty   } from '../../modules/utils/functions.nf'
include {assertKeyExistsAndNotEmpty } from '../../modules/utils/functions.nf'
include {FASTQC                     } from '../../modules/fastqc'
include {FASTQC as FASTQC_RAW       } from '../../modules/fastqc'
include {FASTQC as FASTQC_FILTERED  } from '../../modules/fastqc'
include {MULTIQC                    } from '../../subworkflows/multiqc'
include {META_TO_PED                } from '../../subworkflows/pedigree/meta2ped'
include {FASTQ_HEAD                 } from '../../modules/fastq/head'
include {BWA_INDEX                  } from '../../modules/bwa/index'
include {BWA_MEM                    } from '../../modules/bwa/mem'
include {FASTP                      } from '../../modules/fastp'
include { PREPARE_ONE_REFERENCE     } from '../../subworkflows/samtools/prepare.one.ref'
include { SAMTOOLS_STATS            } from '../../modules/samtools/stats'
include { UNMAPPED                  } from '../../subworkflows/unmapped'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}




workflow {
	if(params.samplesheet==null) {
		log.error("undefined --samplesheet");
		exit -1;
		}
    metadata=["id":"fastqc"];
	versions = Channel.empty()
	multiqc = Channel.empty()
	ch0 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true,sep:",")
        
    META_TO_PED(metadata,ch0)
    versions = versions.mix(META_TO_PED.out.versions)    
        
        
    ch1 = ch0.map{assertKeyExistsAndNotEmpty(it,"sample")}
        .map{
            if(it.fastq_1==null && it.R1!=null) return it.plus(fastq_1:it.R1);
            return it;
            }
        .map{
            if(it.fastq_2==null && it.R2!=null) return it.plus(fastq_2:it.R2);
            return it;
            }  
        .map{assertKeyExistsAndNotEmpty(it,"fastq_1")}
        .branch{
                paired: it.fastq_2!=null &&  !it.fastq_2.isEmpty() &&  !it.fastq_2.equals(".")
                single: true
                }
    
    
    
    
    ch2a = ch1.paired.map{assertKeyExistsAndNotEmpty(it,"fastq_2")}.map{[
        [id:it.sample],
        [file(it.fastq_1),file(it.fastq_2)]
        ]}

    ch2b = ch1.single.map{[
        [id:it.sample],
        [file(it.fastq_1)]
        ]}
	
    /* small bwa */
    if(params.fasta!=null && (params.fastq_head_N as int)>0) {
        PREPARE_ONE_REFERENCE(
            metadata.plus(skip_scatter:true),
            Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
            )
        versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
        if(params.bwa_mem_dir==null) {
            BWA_INDEX(PREPARE_ONE_REFERENCE.out.fasta)
            versions = versions.mix(BWA_INDEX.out.versions)
            bwa_mem_dir =  BWA_INDEX.out.bwa_index
            }
        else
            {
            bwa_mem_dir = [[id:"userbwamemdir"],file(params.bwa_mem_dir)] 
            }

        FASTQ_HEAD(
            ch2a
                .flatMap{meta,R1R2->[
                    [meta.plus(id:meta.id+"_R1"),R1R2[0]],
                    [meta.plus(id:meta.id+"_R2"),R1R2[1]]
                    ]}
                .mix(ch2b)
            )
        versions = versions.mix(FASTQ_HEAD.out.versions)

        FASTQC_RAW(FASTQ_HEAD.out.fastq.map{meta,fq->[meta,[fq]]})
        versions = versions.mix(FASTQC_RAW.out.versions)
        multiqc = multiqc.mix(FASTQC_RAW.out.zip)
        multiqc = multiqc.mix(FASTQC_RAW.out.html)

        FASTP(FASTQ_HEAD.out.fastq.map{meta,fq->[meta,[fq]]})
        versions = versions.mix(FASTP.out.versions)
        multiqc = multiqc.mix(FASTP.out.json)

        FASTQC_FILTERED(FASTP.out.fastqs)
        versions = versions.mix(FASTQC_FILTERED.out.versions)
        multiqc = multiqc.mix(FASTQC_FILTERED.out.zip)
        multiqc = multiqc.mix(FASTQC_FILTERED.out.html)

        BWA_MEM(
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            bwa_mem_dir,
            [[id:"nobed"],[]],
            FASTP.out.fastqs
                .map{meta,fqs->[meta,(fqs instanceof List?fqs:[fqs])]}
                .map{meta,fqs->[meta,fqs[0],[]]}
            )
        versions = versions.mix(BWA_MEM.out.versions)


        SAMTOOLS_STATS(
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            [[id:"nobed"],[]],
            BWA_MEM.out.bam
            )
        versions = versions.mix(SAMTOOLS_STATS.out.versions)
        multiqc = multiqc.mix(SAMTOOLS_STATS.out.stats)

        UNMAPPED(
            metadata,
            [[id:"contaminants"],file("${moduleDir}/../unmapped/contaminants.txt")],
            BWA_MEM.out.bam
		.combine(PREPARE_ONE_REFERENCE.out.fasta)
		.combine(PREPARE_ONE_REFERENCE.out.fai)
		.combine(PREPARE_ONE_REFERENCE.out.dict)
		.map{meta,bam,bai,_m2,fasta,_m3,fai,_m4,dict->[meta,bam,bai,fasta,fai,dict]}
            )
        versions = versions.mix(UNMAPPED.out.versions)
        multiqc = multiqc.mix(UNMAPPED.out.multiqc)
        }
    
	

    FASTQC( ch2a.mix(ch2b) )
    versions = versions.mix(FASTQC.out.versions)
    multiqc = multiqc.mix(FASTQC.out.zip)
    multiqc = multiqc.mix(FASTQC.out.html)

    MULTIQC(
		metadata,
		META_TO_PED.out.sample2collection,
		versions,
        [[id:"mqc_config"],file("${moduleDir}/multiqc_config.yaml")],
		multiqc
		)
    }


runOnComplete(workflow);

