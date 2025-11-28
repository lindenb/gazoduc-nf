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
nextflow.enable.dsl=2


include {ULTRA_RARES_ITERATION as ITER_FIRST       } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_SECOND      } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_THIRD       } from './iteration.part.nf'
include { UNMAPPED                                 } from '../../subworkflows/unmapped'
include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { SAMTOOLS_STATS                           } from '../../modules/samtools/stats'
include { MULTIQC                                  } from '../../modules/multiqc'
include { COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'



if( params.help ) {
    gazoduc.usage().
        name("Ultrarares").
        desc("ultrares").
        print();
    exit 0
} else {
   gazoduc.validate();
}


List fractionate(List L) {
	L = new ArrayList<>(L);
	def L10 = [];
	int i;
	while(i< L.size())
		def item = L[i];
		if(!isBlank(item[0].status) && item[0].status!="case") {
			L10.add(item);
			L.remove(i);
			}
		else
			{
			++i;
			}
		}
	// sort using biggest files first, exomes at the end
	Collections.sort(L,(A,B)->Long.compare(B[1].size(),A[1].size()))
	while(L10.size()<10 && !L.isEmpty()) {
		L10.add(L.remove(0));
		}

	}


workflow {
	def metadata = [id:"ultrares"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
			throw new IllegalArgumentException("undefined --samplesheet");
			}

  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata,
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/* no fastq samplesheet */
	READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        )
	versions = versions.mix(READ_SAMPLESHEET.out.versions)

	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)

	
	ch1 = META_TO_BAMS.out.bams
		.toArray()
		.flatMap{fractionate(it)}
		.branch{v->
			list10 : v[0] 
			list100: v[1]
			list1000: v[2]
			}

	
	ITER_FIRST(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			bed,
			ch1.list10
			)
    version_ch = version_ch.mix(ITER_FIRST.out.versions)

	iter2_ch  = ITER_SECOND([split_vcf_method:" --variants-count  10000 "], genomeId, iter1_ch.vcf, splitbams_ch.output2, bed)
        version_ch = version_ch.mix(iter2_ch.version)

	iter3_ch  = ITER_THIRD([split_vcf_method:" --variants-count  100 "], genomeId, iter2_ch.vcf, splitbams_ch.output3, bed)
        version_ch = version_ch.mix(iter3_ch.version)

        version_ch = MERGE_VERSION("rare-gatk",version_ch.collect())

emit:
        vcf = iter1_ch.vcf
        index = iter1_ch.index
        version= version_ch

}
