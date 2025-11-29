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
nextflow.enable.dsl=2

include {runOnComplete} from '../../../modules/utils/functions.nf'
//include {ANNOTATE_VCF_01} from '../../../subworkflows/annotation/annotation.vcf.01.nf'
//include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
//include {JVARKIT_VCF_TO_INTERVALS_01} from '../../../subworkflows/jvarkit/vcf2intervals/main.nf'
//include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
//include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
//include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed/main.nf'
//include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GTF   } from '../../../modules/gtf/download/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GFF3  } from '../../../modules/gtf/download/main.nf'
//include {ANNOTATE                               } from '../../../subworkflows/annotation/annotsnv/main.nf'
include {ANNOTATE                                } from '../../../subworkflows/annotation/annotsnv/main.nf'
include {SPLIT_VCF                               } from '../../../subworkflows/jvarkit/splitnvariants/main.nf'
include {BCFTOOLS_CONCAT                         } from '../../../modules/bcftools/concat/main.nf'
include { VCF_INPUT as INPUT_VCF                 } from '../../../subworkflows/nf/vcf_input' 
include {PREPARE_ONE_REFERENCE                   } from '../../../subworkflows/samtools/prepare.one.ref'


workflow {
	versions = Channel.empty()
	def metadata = [id:"annot"]


    if(params.fasta==null) {
      throw new IllegalArgumentException("undefined --fasta");
      }
 	if(params.input==null) {
      throw new IllegalArgumentException("undefined --input");
      }



  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	def pedigree = [ ref_hash, (params.pedigree ? file(params.pedigree) : []) ]
	def bed =      [ ref_hash, (params.bed ? file(params.bed) : []) ]
	
	
  /***************************************************
   *
   *  INPUT_VCF
   *
   */
 	INPUT_VCF(
      metadata.plus([
        arg_name : "input",
        path: params.input,
        require_index : true,
        required : true,
        unique : true
        ])
      )
    versions = versions.mix(INPUT_VCF.out.versions)
    vcfs = INPUT_VCF.out.vcf

	

	 /** break the VCF into parts */
    SPLIT_VCF(
        [id:"annot"],
        fasta,
        fai,
        dict,
        bed,
        vcfs
        )
    vcfs = SPLIT_VCF.output.vcf
	
	DOWNLOAD_GTF(dict)
	DOWNLOAD_GFF3(dict)
	
	ANNOTATE(
		[id:"annot"],
		fasta,
		fai,
		dict,
		DOWNLOAD_GTF.out.output,
		DOWNLOAD_GFF3.out.output,		
		pedigree,
		vcfs
		)

	BCFTOOLS_CONCAT(
		ANNOTATE.out.vcf
			.map{[ it[0], [it[1],it[2]] ]}
			.groupTuple()
			.map{[it[0],it[1].flatten()]} , 
		[[id:"nobed"],[]]
		)

	versions= versions.mix(ANNOTATE.out.versions)
	}


runOnComplete(workflow);
