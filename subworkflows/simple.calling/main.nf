/*

Copyright (c) 2024 Pierre Lindenbaum

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
include {HAPLOTYPECALLER_DIRECT  } from '../../modules/gatk/hapcallerdirect'
include {BCFTOOLS_CALL           } from '../../modules/bcftools/call'
include {BCFTOOLS_MERGE          } from '../../modules/bcftools/merge3'
include {BCFTOOLS_CONCAT         } from '../../modules/bcftools/concat3'
include {CALL_FREEBAYES          } from '../../modules/freebayes/call'
include {isBlank                 } from '../../modules/utils/functions'

/** universal calling */
workflow SIMPLE_CALLING {
	take:
		metadata
		fasta
		fai
		dict
        dbsnp
		pedigree
		beds
		bams //[ meta, bam,bai]
	main:
		versions = Channel.empty()


        if(metadata.method==null) {
            log.warn("CALLING metadata.method is mssing");
			throw new IllegalArgumentException("In SIMPLE_CALLING bed should have a unique id");
            }
         if(!metadata.method.matches("(gatk|freebayes|bcftools)")) {
            log.warn("CALLING metadata.method is unknown: ${metadata.method}");
			throw new IllegalArgumentException("In SIMPLE_CALLING unknwon method ${metadata.method}");
            }

		/* checl all beds have an unique ID */
		beds.map{meta,bed->[meta.id,bed]}
			.groupTuple()
			.filter{it[1].size()!=1}
			.map{throw new IllegalArgumentException("In CALLING: bed should have a unique id: ${it}"); return it;}
	

		/* group all the bams, make a (meta,bams,bed] */
		ch2 = bams.map{meta,bam,bai->["${meta.batch?:meta.id}", [bam, bai]]}
			.groupTuple()
			.map{group_id,bam_files->[[batch:group_id],bam_files.flatten().sort()]}
			.combine(beds)
			.map{meta1,bam_files,meta2,bed->[meta1.plus(id:meta2.id),bam_files,bed]} //give them the bed.id

        vcf_out= Channel.empty()

        if(metadata.method=="bcftools") {
			BCFTOOLS_CALL(
				fasta,
				fai,
				pedigree,
				[[id:"noploidy"],[]],
				ch2
				)
			versions = versions.mix(BCFTOOLS_CALL.out.versions)
			vcf_out = BCFTOOLS_CALL.out.vcf
			}
        else if(metadata.method=="freebayes") {
            CALL_FREEBAYES(
                    fasta,
                    fai,
                    ch2
                    ) 
			versions = versions.mix(CALL_FREEBAYES.out.versions)
			vcf_out = CALL_FREEBAYES.out.vcf
            }
		else if(metadata.method=="gatk") {
			HAPCALLER(
				fasta,
				fai,
				dict,
				dbsnp,
				pedigree,
				ch1
				)
			versions = versions.mix(metadata.out.versions)
			vcf_out = metadata.out.vcf
            }
        else {
            throw new IllegalArgumentException("In CALLING: unknown method ${metadata.method}");
            }

		/* group data by BATCH */
		group_by_batch = vcf_out
				.map{meta,vcf,tbi,bed->[meta.batch,[vcf,tbi]]} // group by bedid
				.groupTuple()
				.map{batch_id,vcf_files->[[id:batch_id],vcf_files.flatten().sort()]}
				
		/** concat per batch, so we can have one vcf per batch/sample  */
		BCFTOOLS_CONCAT(group_by_batch)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)
		per_group_vcf = BCFTOOLS_CONCAT.out.vcf

		if(metadata.with_merge==null || parseBoolean(metadata.with_merge)) {
			// merge those having more than one sample
			BCFTOOLS_MERGE(
				BCFTOOLS_CONCAT.out.vcf
					.map{meta,vcf,tbi-> [ [id:"bcftools.call"],[vcf,tbi]]}
					.groupTuple()
					.map{meta,files ->[meta, files.flatten().sort() ]}
				)
			versions = versions.mix(BCFTOOLS_MERGE.out.versions)
			vcf_out = BCFTOOLS_MERGE.out.vcf
			}
		else
			{
			vcf_out = Channel.empty()
			}
		
	emit:
		versions
		per_group_vcf
		vcf = vcf_out
	}


