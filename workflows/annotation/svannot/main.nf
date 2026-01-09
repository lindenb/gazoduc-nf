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


include { validateParameters                          } from 'plugin/nf-schema'
include { paramsHelp                                  } from 'plugin/nf-schema'
include { paramsSummaryLog                            } from 'plugin/nf-schema'
include { samplesheetToList                           } from 'plugin/nf-schema'
include { runOnComplete;dumpParams                    } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                       } from '../../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                                   } from '../../../subworkflows/nf/vcf_input'
include { DOWNLOAD_GTF_OR_GFF3  as  DOWNLOAD_GTF      } from '../../../modules/gtf/download'
include { DOWNLOAD_GNOMAD_SV                          } from '../../../modules/gnomad_sv/download.vcf'
include { GTF_TO_BED                                  } from '../../../modules/jvarkit/gtf2bed'
include { JVARKIT_VCFGNOMADSV                         } from '../../../modules/jvarkit/vcfgnomadsv'
include { BEDTOOLS_SLOP                               } from '../../../modules/bedtools/slop'
include { BEDTOOLS_MERGE                              } from '../../../modules/bedtools/merge'
include { BCFTOOLS_VIEW                               } from '../../../modules/bcftools/view'
include { ANNOTSV                                     } from '../../../subworkflows/annotsv'


workflow {
		validateParameters()

		if( params.help ) {
			log.info(paramsHelp())
			exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}
	multiqc = Channel.empty()
	versions = Channel.empty()
	def metadata = [id:"annotsv"]
		

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	
	DOWNLOAD_GTF(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(DOWNLOAD_GTF.out.versions)

    DOWNLOAD_GNOMAD_SV( PREPARE_ONE_REFERENCE.out.dict )
    versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)


	VCF_INPUT(
        metadata.plus([
			path:params.vcf,
			require_index:true,
			arg_name:"vcf",
			required: true,
			unique : false
			])
        )
    versions = versions.mix(VCF_INPUT.out.versions)

	if(params.bed==null) {
		bed = PREPARE_ONE_REFERENCE.out.scatter_bed
	} else if(params.bed.matches("(coding[_-])?(exon|gene)s?")) {
		GTF_TO_BED(
			PREPARE_ONE_REFERENCE.out.dict,
			DOWNLOAD_GTF.out.gtf.map{meta,gtf,_tbi->[meta,gtf]}
			)
		versions = versions.mix(GTF_TO_BED.out.versions)
		bed =  GTF_TO_BED.out.bed
	} else {
		bed =Channel.of(params.bed).map{[[id:"bed"],file(f)]}
	}

	BEDTOOLS_SLOP(
		PREPARE_ONE_REFERENCE.out.fai,
		bed
		)
	versions = versions.mix(BEDTOOLS_SLOP.out.versions)

	BEDTOOLS_MERGE( BEDTOOLS_SLOP.out.bed )
	versions = versions.mix(BEDTOOLS_MERGE.out.versions)

	BCFTOOLS_VIEW(
		BEDTOOLS_MERGE.out.bed,
		[[id:"no_sample"],[]],
		VCF_INPUT.out.vcf
		)
	versions = versions.mix(BCFTOOLS_VIEW.out.versions)

	JVARKIT_VCFGNOMADSV(
		DOWNLOAD_GNOMAD_SV.out.vcf,
		BCFTOOLS_VIEW.out.vcf
		)
	versions = versions.mix(JVARKIT_VCFGNOMADSV.out.versions)


	ANNOTSV(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		JVARKIT_VCFGNOMADSV.out.vcf
		)
	versions = versions.mix(ANNOTSV.out.versions)
}