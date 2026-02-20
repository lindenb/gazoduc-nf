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

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
  // log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
//validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)

include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { VQSR                                     } from '../../subworkflows/vqsr'
include { MULTIQC                                  } from '../../modules/multiqc'
include { COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include { VCF_STATS                                } from '../../subworkflows/vcfstats'
include { runOnComplete; dumpParams                } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'



if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
		versions = Channel.empty()
		multiqc = Channel.empty()

		def metadata = [  id: "vqsr"  ]
		
		if(params.fasta==null) {
			log.error("--fasta undefined")
			exit -1;
			}
		if(params.vcf==null) {
			log.error("--vcf undefined")
			exit -1;
			}
		if(params.recal_snps_args==null) {
			log.error("--recal_snps_args undefined")
			exit -1;
			}
		if(params.recal_indels_args==null) {
			log.error("--recal_indels_args undefined")
			exit -1;
			}

		PREPARE_ONE_REFERENCE(
			metadata.plus(skip_scatter:true),
			Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
			)
		versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

		def dbsnp  = [[id:"npdbsnp"], [] , [] ] 
		if(params.dbsnp!=null) {
			dbsnp  = [[id:"dbsnp"], file(params.dbsnp) , file(params.dbsnp+".tbi") ] 
			}
		
		def recal_snps = [[id:"snps"], params.recal_snps_args]
		def recal_indels = [[id:"indels"],params.recal_indels_args]

		VCF_INPUT(metadata.plus([
			path: params.vcf,
			arg_name: "vcf",
			require_index : true,
			required: true,
			unique : false
			]))
		versions = versions.mix(VCF_INPUT.out.versions)

		
	
		VQSR(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			dbsnp,
			recal_snps,
			recal_indels,
			in_vcf_ch
			)
		versions = versions.mix(VQSR.out.versions)

	if(params.with_stats==true) {

		def gtf    = [metadata, file(params.gtf), file(params.gtf+".tbi") ]
		def gff3    = [metadata, file(params.gff3), file(params.gff3+".tbi") ]


		VCF_STATS(
				metadata,
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
				Channel.of(gtf),
				Channel.of(gff3),
				[[id:"noped"],[]],
				[[id:"nogroup2sample"],[]],
				[[id:"nobed"],[]],
				VQSR.out.vcf
				)
		versions = versions.mix(VCF_STATS.out.versions)
		multiqc = versions.mix(VCF_STATS.out.multiqc)
		}

	if(params.with_multiqc==true) {
		COMPILE_VERSIONS(versions.collect())
		multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

		MULTIQC(multiqc.map{meta,f->f}.collect().map{[[id:"vqsr"],it]})
		}
	}
runOnComplete(workflow)
