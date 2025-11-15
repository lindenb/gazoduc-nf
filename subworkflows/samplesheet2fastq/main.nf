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
include { isBlank                    } from '../../modules/utils/functions'
include { assertKeyExistsAndNotEmpty } from '../../modules/utils/functions'
include { ORA_TO_FASTQ               } from '../../subworkflows/ora/ora2fastq'
include { BAM_TO_FASTQ               } from '../../subworkflows/bam2fastq'
include { SAMTOOLS_SAMPLES           } from '../../modules/samtools/samples'
include { isEmptyGz                  } from '../../modules/utils/functions.nf'
include { extractIlluminaName        } from '../../modules/utils/functions.nf'
include { verify                     } from '../../modules/utils/functions.nf'
include { parseBoolean               } from '../../modules/utils/functions.nf'
include { FASTQ_MERGE as MERGE_SE    } from '../../modules/fastq/merge'
include { FASTQ_MERGE as MERGE_PE1   } from '../../modules/fastq/merge'
include { FASTQ_MERGE as MERGE_PE2   } from '../../modules/fastq/merge'

boolean hasKey(def h, def id) {
	return h!=null && !isBlank(h[id]);
	}
Map cleanupHash(Map h) {
	return h.findAll{k,v->!k.matches("fasta|fai|dict|bam|bai|ora|R1|R2|fastq_1|fastq_2|bed")}
	}



workflow SAMPLESHEET_TO_FASTQ {
take:
    workflow_metadata // bam2fastq_method=(samtools|gatk|parabricks)
    samplesheet // row from splitCsv...
main:
    versions = Channel.empty()
    paired = Channel.empty()
    single  = Channel.empty()
    
    if(workflow_metadata.merge_fastqs==null) {
        log.warn("SAMPLESHEET_TO_FASTQ: merge_fastqs undefined.")
        }


    ch0 = samplesheet
        .branch {
        	fastq :  hasKey(it,"fastq_1") && !hasKey(it,"bam") && !it.fastq_1.endsWith(".ora")
        	ora   :  hasKey(it,"fastq_1") && !hasKey(it,"bam") &&  it.fastq_1.endsWith(".ora") 
        	bam   : !hasKey(it,"fastq_1") &&  hasKey(it,"bam")
        	other : true
        	}
     
    ch0.other.map{throw new IllegalArgumentException("SAMPLESHEET_TO_FASTQ: undefined/ambigous input ${it}.");}
    
    /******************************************
     *
     * INPUT is ORA
     *
     */
    ora_ch =  ch0.ora.map{
            verify( it.fastq_1.endsWith(".ora") , "R1 should end with ORA ${it}" );
            verify( isBlank(it.fastq_2) || it.fastq_2.endsWith(".ora") , "SAMPLESHEET_TO_FASTQ: both R1 and R2 should end with '.ora':  ${it}.");
            return it;
            }
        .map{
            if(hasKey(it,"id")) return it;
            if(hasKey(it,"sample")) return it.plus(id:it.sample);
            def ilm_name_1 = extractIlluminaName(it.fastq_1);
            if(ilm_name_1!=null) {
                //verify( ilm_name_1.side=="R1" , "SAMPLESHEET_TO_FASTQ: Expected side:R1 :  ${it} ."); TODO fix rotavirus path
                // if paired end check R2 share the same info with R1
                if(!isBlank(it.fastq_2)) {
                    def ilm_name_2 = extractIlluminaName(it.fastq_2);
                    verify( ilm_name_2!=null , "SAMPLESHEET_TO_FASTQ: can extract ilmn name for R1 but not from R2 ?? :  ${it}.");
                    verify( ilm_name_1.id == ilm_name_2.id , "SAMPLESHEET_TO_FASTQ: Discordant illumina names :  ${it}.");
                    verify( ilm_name_1.lane == ilm_name_2.lane , "SAMPLESHEET_TO_FASTQ: Discordant illumina lanes :  ${it}.");
                    verify( ilm_name_1.index == ilm_name_2.index , "SAMPLESHEET_TO_FASTQ: Discordant illumina index :  ${it}.");
                    verify( ilm_name_1.split == ilm_name_2.split , "SAMPLESHEET_TO_FASTQ: Discordant illumina index :  ${it}.");
                    verify( ilm_name_2.side == "R2" ,  "SAMPLESHEET_TO_FASTQ: Expected side:R2 :  ${it}.");
                    }
                return it.plus(id:ilm_name_1.id);
                }
            verify(!isBlank(it.id),"SAMPLESHEET_TO_FASTQ: id missing for  ${it}.");
            }
        .map{
            def hash =it.findAll{k,v->!k.matches("(fastq_1|fastq_2)")}
            if(isBlank(it.fastq_2)) {
                return [hash,file(it.fastq_1)]
                }
            return [hash , file(it.fastq_1), file(it.fastq_2) ];
            }
    	
    ORA_TO_FASTQ(
     	workflow_metadata,
     	ora_ch
     	)
    versions = versions.mix(ORA_TO_FASTQ.out.versions)
    paired = paired.mix(ORA_TO_FASTQ.out.paired_end)
    single = single.mix(ORA_TO_FASTQ.out.single_end)





    
    /******************************************
     *
     * INPUT is BAM OR CRAM
     *
     */
    ch3 = ch0.bam.branch {v->
        no_sample:  !hasKey(v,"sample") && !hasKey(v,"id")
        has_sample : true
        }
    
    /* collect all references */
    all_references = ch0.bam
        .filter{hasKey(it,"fasta")}
        .map{
            if(hasKey(it,"fai")) return it;
            return it.plus(fai: it.fasta+".fai")
            }
        .flatMap{[it.fasta,it.fai]}
        .unique()
        .map{file(it)}
        .toSortedList() // enable empty array , not 'collect()'
        .map{[[id:"references"],it]}

    //undefined Sample name for bam, call samtools samples to resole sample name
    SAMTOOLS_SAMPLES(
        all_references,
        ch3.no_sample
            .map{file(it.bam)}
            .collect()
            .map{[workflow_metadata,it.sort()]}
        )
    versions = versions.mix(SAMTOOLS_SAMPLES.out.versions)
    

    resolved_bam_sample_ch = SAMTOOLS_SAMPLES.out.samplesheet
        .map{it[1]}
        .splitCsv(header:false,sep:'\t')
        .map{it->[it[0],it[1],it[2]]} //extract sample bam  and ref
        .map{sample,bam,ref->
            if(isBlank(sample)) throw new IllegalArgumentException("no sample/id defined for ${bam}");
            return [sample,bam,ref];
            }
        .map{sample,bam,ref->[file(bam).toRealPath().toString(),sample,ref]}
        .join(ch3.no_sample.map{meta->[file(meta.bam).toRealPath().toString(),meta]}, failOnMismatch:true)
        .map{bam,sn,ref,meta->
            def h = meta.plus([id:sn,sample:sn])
            if(isBlank(h.fasta) && !isBlank(ref)) {
                h = h.plus([fasta:ref, fai: ref+".fai"]);
                }
            return h;
            }

    

     BAM_TO_FASTQ(
        workflow_metadata,
     	ch3.has_sample
            .mix(resolved_bam_sample_ch)
            .map{                
                if(!hasKey(it,"fasta") && it.bam.endsWith(".cram"))  throw new IllegalArgumentException("no fasta sequence defined for ${it.bam}"); 
                return it;
                }
            .map{
                if(hasKey(it,"id")) return it;
                if(hasKey(it,"sample")) return it.plus(id:it.sample);
                throw new IllegalArgumentException("undefined id in ${it}");
                }
            .map{assertKeyExistsAndNotEmpty(it,"id")}
            .map{[
                cleanupHash(it),
                file(it.bam),
                (hasKey(it,"bai")?file(it.bai): file(it.bam+(it.bam.endsWith(".bam")?".bai":".crai"))),
                (hasKey(it,"fasta")?file(it.fasta):[]),
                (hasKey(it,"fai")?file(it.fai): (hasKey(it,"fasta")?file(it.fasta+".fai"):[])),
                (hasKey(it,"dict")?file(it.dict):(hasKey(it,"fasta")?file(it.fasta.replaceAll("\\.(fasta|fa|fna)\$",".dict")):[])),
                (hasKey(it,"bed")?file(it.bed): [] )
                ]}
     	)
    versions = versions.mix(BAM_TO_FASTQ.out.versions)
    paired = paired.mix(BAM_TO_FASTQ.out.paired_end) 
    single = single.mix(BAM_TO_FASTQ.out.single_end) 


     /**
      * REGULAR FASTQ
      */
    
    ch1 =  ch0.fastq
        .map{
            if(hasKey(it,"id")) return it;
            if(hasKey(it,"sample")) return it.plus(id:it.sample);
            verify(!isBlank(it.fastq_1),"fastq_1 missing in ${it}");

            def ilm_name_1 = extractIlluminaName(it.fastq_1);
            if(ilm_name_1!=null) {
                verify( ilm_name_1.side=="R1" , "SAMPLESHEET_TO_FASTQ: Expected side:R1 :  ${it}.");
                // if paired end check R2 share the same info with R1
                if(!isBlank(it.fastq_2)) {
                    def ilm_name_2 = extractIlluminaName(it.fastq_2);
                    verify( ilm_name_2!=null , "SAMPLESHEET_TO_FASTQ: can extract ilmn name for R1 but not from R2 ?? :  ${it}.");
                    verify( ilm_name_1.id == ilm_name_2.id , "SAMPLESHEET_TO_FASTQ: Discordant illumina names :  ${it}.");
                    verify( ilm_name_1.lane == ilm_name_2.lane , "SAMPLESHEET_TO_FASTQ: Discordant illumina lanes :  ${it}.");
                    verify( ilm_name_1.index == ilm_name_2.index , "SAMPLESHEET_TO_FASTQ: Discordant illumina index :  ${it}.");
                    verify( ilm_name_1.split == ilm_name_2.split , "SAMPLESHEET_TO_FASTQ: Discordant illumina index :  ${it}.");
                    verify( ilm_name_2.side == "R2" ,  "SAMPLESHEET_TO_FASTQ: Expected side:R2 :  ${it}.");
                    }
                return it.plus(id:ilm_name_1.id);
                }
            verify(!isBlank(it.id),"id missing in ${it}");
            return it;
            }
        .map{assertKeyExistsAndNotEmpty(it,"id")} 
        .map{assertKeyExistsAndNotEmpty(it,"fastq_1")}
        .branch{
                paired: !isBlank(it.fastq_1) && !isBlank(it.fastq_2)
                single: true
                }

        paired = paired.mix(
            ch1.paired
                .map{
                    if(it.fastq_1.equals(it.fastq_2)) {
                        throw new IllegalArgumentException("R1==R2 in ${it}.");
                        }
                    return it;
                    }
                .map{[
                    cleanupHash(it),
                    file(it.fastq_1),
                    file(it.fastq_2)
                    ]}
            )
        
        single = single.mix(ch1.single
                .map{[
                    cleanupHash(it),
                    file(it.fastq_1)
                    ]}
                )
        /* file might not exist in stub mode */
        if(!workflow.stubRun) {
            paired = paired.filter{meta,R1,R2->!(isEmptyGz(R1) && isEmptyGz(R2))}
            single = single.filter{meta,R1->!(isEmptyGz(R1))}
            }

        if(workflow_metadata.merge_fastqs!=null && parseBoolean(workflow_metadata.merge_fastqs)) {
            /** merge single end */
            ch4 = single
                .map{meta,R1->[meta.id,meta,R1]}
                .groupTuple()
                .map{_meta_id,metas,files->[metas[0],files.sort()]}
                .branch{meta,files->
                    one: files.size()==1
                    multi: true
                    }
           
            MERGE_SE(ch4.multi.map{meta,files->[meta.plus(prefix:meta.id+".R0"),files]})
            versions = versions.mix(MERGE_SE.out.versions)

            single = ch4.one
                    .map{meta,files->[meta,files[0]]}
                    .mix(MERGE_SE.out.fastq)
            
            /** merge paired-end */
            ch5 = paired
                .map{meta,R1,R2->[meta.id,meta,R1,R2]}
                .groupTuple()
                .map{_meta_id,metas,files1,files2->[metas[0],files1.sort(),files2.sort()]}
                .branch{meta,files1,files2->
                    one: files1.size()==1
                    multi: true
                    }
            MERGE_PE1(ch5.multi.map{meta,f1,f2->[meta.plus(prefix:meta.id+".R1"),f1]})
            versions = versions.mix(MERGE_PE1.out.versions)

            MERGE_PE2(ch5.multi.map{meta,f1,f2->[meta.plus(prefix:meta.id+".R2"),f2]})
            versions = versions.mix(MERGE_PE2.out.versions)
        

            paired = ch5.one
                    .map{meta,f1,f2->[meta,f1[0],f2[0]]}
                    .mix(
                        MERGE_PE1.out.fastq
                            .map{meta,f->[meta.findAll{k,v->!k.matches("prefix")},f]}
                            .join(MERGE_PE2.out.fastq.map{meta,f->[meta.findAll{k,v->!k.matches("prefix")},f]})
                        )
            }

		paired.mix(single).filter{
			verify(it[0].id.matches("[A-Za-z0-9][A-Za-z_0-9\\.\\-]*"),"Bad id for a fastq ${it}");
			}


    emit:
        versions
        paired_end = paired
        single_end = single
}
