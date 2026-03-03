
/*

THIS FILE WAS GENERATED DO NOT EDIT !!!

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

include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { APHAGENOME_VARIANT_SCORING               } from '../../subworkflows/alphagenome/variant_scoring'

workflow {
  versions = Channel.empty()
  multiqc = Channel.empty()


	def metadata = [id:"alphagenomevcf"]
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

  /***************************************************
   *
   *  LOAD VCF
   *
   */
	VCF_INPUT(metadata.plus([
		path: params.vcf,
		arg_name: "vcf",
		require_index : true,
		required: true,
		unique : false
		]))
	versions = versions.mix(VCF_INPUT.out.versions)

    APHAGENOME_VARIANT_SCORING(
        metadata,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.dict,
        PREPARE_ONE_REFERENCE.out.fai,
        VCF_INPUT.out.vcf
        )
	versions = versions.mix(APHAGENOME_VARIANT_SCORING.out.versions)


}

