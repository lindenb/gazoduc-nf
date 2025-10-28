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
    
    ch0 = samplesheet
        .branch {
        	fastq :  hasKey(it,"fastq_1") && !hasKey(it,"ora") && !hasKey(it,"bam")
        	ora   : !hasKey(it,"fastq_1") &&  hasKey(it,"ora") && !hasKey(it,"bam")
        	bam   : !hasKey(it,"fastq_1") && !hasKey(it,"ora") &&  hasKey(it,"bam")
        	other : true
        	}
     
    ch0.other.map{throw new IllegalArgumentException("SAMPLESHEET_TO_FASTQ: undefined/ambigous input ${it}.");}
    
    /******************************************
     *
     * INPUT is ORA
     *
     */
     
    ORA_TO_FASTQ(
     	workflow_metadata,
     	ch0.ora
            .map{
                if(hasKey(it,"id")) return it;
                if(hasKey(it,"sample")) return it.plus(id:it.sample);
                throw new IllegalArgumentException("undefined id in ${it}");
                }
            .map{assertNotEmpty(it,"id")}
            .map{[
                cleanupHash(it),
                file(it.ora)
                ]}
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
            if(isBlank(ref) && bam.endsWith(".cram"))  throw new IllegalArgumentException("no fasta sequence defined for ${bam}");
            return it;
            }
        .map{sample,bam,ref->[file(bam).toRealPath().toString(),sample,ref]}
        .join(ch3.no_sample.map{meta->[file(meta.bam).toRealPath().toString(),meta]}, failOnMismatch:true)
        .map{bam,sn,ref,meta->
            def h = meta.plus([id:sample,sample:sn])
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
                if(hasKey(it,"id")) return it;
                if(hasKey(it,"sample")) return it.plus(id:it.sample);
                if(!hasKey(it,"fasta") && it.bam.endsWith(".cram"))  throw new IllegalArgumentException("no fasta sequence defined for ${it.bam}");
                throw new IllegalArgumentException("undefined id in ${it}");
                }
            .map{assertNotEmpty(it,"id")}
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
    
    ch1 = ch0.fastq
        .map{
            if(hasKey(it,"id")) return it;
            if(hasKey(it,"sample")) return it.plus(id:it.sample);
            throw new IllegalArgumentException("undefined id in ${it}");
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


    paired = paired.filter{meta,R1,R2->!(isEmptyGz(R1) && isEmptyGz(R2))}
    single = single.filter{meta,R1->!(isEmptyGz(R1))}

emit:
    versions
    paired_end = paired
    single_end = single
}
