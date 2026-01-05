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
include {INDEXCOV_REBUILD_BAI                } from '../../modules/indexcov/rebuild.bai'
include {APPLY_INDEXCOV                      } from '../../modules/indexcov/apply'
include {MERGE_BEDS                          } from '../../modules/indexcov/merge.beds'
include {JVARKIT_INDEXCOV2VCF                } from '../../modules/jvarkit/indexcov2vcf'

workflow INDEXCOV {
     take:
		meta
		fasta
		fai
		dict
		scatter_N_bed
		pedigree
		bams
     main:
        versions = Channel.empty()
		multiqc = Channel.empty()

		
		ch1 = bams.branch{
			ok: it.size()==2
			other: true
			}

		ch1.other.map{throw new IllegalArgumentException("INDEXCOV:expected: meta,bam,ref,fasta,fai got ${it}.")}


		INDEXCOV_REBUILD_BAI(
			fasta,
			fai,
			ch1.ok
			)
		versions = versions.mix(INDEXCOV_REBUILD_BAI.out.versions)

	
		
		if(meta.batch_size==null) {
			log.warn("INDEXCOV: meta.batch_size undefined");
			}
		def collate_size = ((meta.batch_size?:1000) as int)

		bams2_ch = INDEXCOV_REBUILD_BAI.out.bam
			.map{meta,bam,bai->[bam,bai]} //bam,bai
			.collect(flat:false)
			.flatMap{
				def L1=[];
				def L2=[];
				def L3 = it.sort({a,b->a[0].name.compareTo(b[0].name)});
				int i=0;
				while(i < L3.size()) {
					if(L2.size()==collate_size*2 /* because bam and bai */) {
						L1.add(L2);
						L2=[];
						}
					L2.add(L3[i][0]);//bam
					L2.add(L3[i][1]);//bai
					i++;
					}
				if(L2.size()>0) L1.add(L2);
				return L1;
				}
			.map{T->[[id:"batch."+T.collect{V->V.toString()}.join(" ").md5().substring(0,7)], T]}
		
    APPLY_INDEXCOV(fasta,fai,bams2_ch)
	versions = versions.mix( APPLY_INDEXCOV.out.versions)

	MERGE_BEDS(
		APPLY_INDEXCOV.out.bed.
			map{T->T[1]}.
			collect().
			map{[[id:meta.id],it]}
		)
	versions = versions.mix(MERGE_BEDS.out.versions)

	bed_ch = MERGE_BEDS.out.bed



	JVARKIT_INDEXCOV2VCF(
		fasta,fai,dict,
		pedigree,
		bed_ch.map{meta,bed,tbi->[meta,bed]}
		)
	versions = versions.mix(JVARKIT_INDEXCOV2VCF.out.versions)

	emit:
		bed  = bed_ch
		zip = APPLY_INDEXCOV.out.zip
		vcf  = JVARKIT_INDEXCOV2VCF.out.vcf
		multiqc
		versions
	}
