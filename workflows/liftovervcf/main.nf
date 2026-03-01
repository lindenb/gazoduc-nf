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
include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { runOnComplete                            } from '../../modules/utils/functions.nf'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { verify                                   } from '../../modules/utils/functions.nf'
include { LIFTOVER_VCF                             } from '../../subworkflows/liftovervcf'


workflow {

		if(!workflow.stubRun) {
		validateParameters()
		}

	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	def metadata = [id:"liftover"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.vcf==null) {
			throw new IllegalArgumentException("undefined --vcf");
			}


  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata.plus(skip_scatter:true),
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/**
	* prepare VCF inputs :VCF may be composed of multiple VCFs (one per contig)
	*/
	VCF_INPUT(metadata.plus([
		path: params.vcf,
		arg_name: "vcf",
		require_index : true,
		required: true,
		unique : false
		]))
	versions = versions.mix(VCF_INPUT.out.versions)

	LIFTOVER_VCF(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		VCF_INPUT.out.vcf
		)
	versions = versions.mix(LIFTOVER_VCF.out.versions)
	multiqc = multiqc.mix(LIFTOVER_VCF.out.multiqc)
	}

runOnComplete(workflow)
