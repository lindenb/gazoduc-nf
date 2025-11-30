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
include {INDEXCOV_REBUILD_BAI                } from '../../../modules/indexcov/rebuild.bai'
include {APPLY_INDEXCOV                      } from '../../../modules/indexcov/apply'
include {MERGE_BEDS                          } from '../../../modules/indexcov/merge.beds'
include {JVARKIT_INDEXCOV2VCF                } from '../../../modules/jvarkit/indexcov2vcf'

workflow INDEXCOV {
     take:
	 	meta
        fasta
		fai
		dict
		pedigree
		bams /* 		meta,bam,(bai)*/

     main:
        versions = Channel.empty()
		multiqc = Channel.empty()

		ch1 = bams.branch{
			rebuild: it.size()==2 || it[1].name.endsWith(".cram") || (it[0].force_rebuild && it[0].force_rebuild==true)
			no_need_to_rebuild: true
			}

		INDEXCOV_REBUILD_BAI(
			fasta,
			fai,
			ch1.rebuild.map{[it[0],it[1]]} //remove bai, if any
			)
		versions = versions.mix(INDEXCOV_REBUILD_BAI.out.versions)

		def collate_size = ((meta.batch_size?:1000) as int)

		bams2_ch = INDEXCOV_REBUILD_BAI.out.bam.mix(ch1.no_need_to_rebuild)
			.map{[it[1],it[2]]} //bam,bai
			.collect(flat:false)
			.flatMap{
				def L1=[];
				def L2=[];
				def L3 = it.sort({a,b->a[0].name.compareTo(b[0].name)});
				int i=0;
				while(i < L3.size()) {
					if(L2.size()==collate_size*2 /* because bam and bai */) {
						L1.add(L2);
						L2=[]
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
			map{[meta,it]}
		)
	versions = versions.mix(MERGE_BEDS.out.versions)

	JVARKIT_INDEXCOV2VCF(
		fasta,fai,dict,
		pedigree,
		MERGE_BEDS.out.bed.map{[it[0],it[1]]}
		)
	versions = versions.mix(JVARKIT_INDEXCOV2VCF.out.versions)


	emit:
		bed  = MERGE_BEDS.out.bed
		zip = APPLY_INDEXCOV.out.zip
		vcf  = JVARKIT_INDEXCOV2VCF.out.vcf
		versions
		multiqc
	}
