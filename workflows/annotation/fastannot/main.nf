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
include { BCFTOOLS_CONCAT                          } from '../../../modules/bcftools/concat3'
include {runOnComplete                             } from '../../../modules/utils/functions.nf'
include {makeKey                                   } from '../../../modules/utils/functions.nf'
include {flatMapByIndex                            } from '../../../modules/utils/functions.nf'
include {parseBoolean                              } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                                } from '../../../subworkflows/nf/vcf_input'
include { VCF_TO_CONTIGS                           } from '../../../subworkflows/bcftools/vcf2contigs'
include { SPLIT_N_VARIANTS                         } from '../../../modules/jvarkit/splitnvariants'
include { SNPEFF_DOWNLOAD                          } from '../../../modules/snpeff/download'
include { SNPEFF_APPLY                             } from '../../../modules/snpeff/apply'
include { JVARKIT_VCF_FLATTEN                      } from '../../../modules/jvarkit/vcfflatten'
include { JVARKIT_VCFGNOMAD                        } from '../../../modules/jvarkit/vcfgnomad'
include { JVARKIT_VCFSTATS                         } from '../../../modules/jvarkit/vcfstats'
include { GROUP_BY_GENE                            } from '../../../modules/jvarkit/groupbygene'
include { BCFTOOLS_INDEX                           } from '../../../modules/bcftools/index'
include { BCFTOOLS_CONTRAST as CONTRAST1           } from '../../../modules/bcftools/contrast'
include { BCFTOOLS_CONTRAST as CONTRAST2           } from '../../../modules/bcftools/contrast'
include { ALPHAMISSENSE_DOWNLOAD                   } from '../../../modules/alphamissense/download'
include { ALPHAMISSENSE_ANNOTATE                   } from '../../../modules/alphamissense/annotate'
include { COMPILE_VERSIONS                         } from '../../../modules/versions/main.nf'
include { MULTIQC                                  } from '../../../modules/multiqc'


workflow  {
		versions = Channel.empty()
		multiqc = Channel.empty()
		def metadata = [id:"genotyping"]

		PREPARE_ONE_REFERENCE(
			metadata.plus(skip_scatter:true),
			Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
			)
		versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
		
		if(params.vcf.endsWith(".list")) {
			vcf_ch = Channel.fromPath(params.vcf).splitText().map{it.trim()}
			}
		else
			{
			vcf_ch = Channel.of(params.vcf)
			}
		
		
		VCF_INPUT(metadata.plus([
				path: params.vcf,
				arg_name: "vcf",
				require_index : true,
				required: true,
				unique : false
				]))
		versions = versions.mix(VCF_INPUT.out.versions)

		VCF_TO_CONTIGS( metadata,VCF_INPUT.out.vcf )
		versions = versions.mix(VCF_TO_CONTIGS.out.versions)

	
		SPLIT_N_VARIANTS(
			[[id:"no_bed"],[]],
			VCF_TO_CONTIGS.out.vcf
			)
		versions = versions.mix(SPLIT_N_VARIANTS.out.versions)
		
		//vcf2inter_ch.splitCsv(header:false,sep:'\t').view()
		SNPEFF_DOWNLOAD(
				PREPARE_ONE_REFERENCE.out.fai
                        .map{meta,_fai->[meta,file(params.snpeff_database_directory)]}
				)
		versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

		SNPEFF_APPLY(
			SNPEFF_DOWNLOAD.out.database,
			SPLIT_N_VARIANTS.out.vcf.flatMap{flatMapByIndex(it,1)}.map{_meta,vcf->[[id:makeKey(vcf)],vcf]}
			)
		versions = versions.mix(SNPEFF_APPLY.out.versions)
		vcf = SNPEFF_APPLY.out.vcf

		JVARKIT_VCFGNOMAD(
			[[id:"gnomad"],file(params.gnomad),file(params.gnomad+".tbi")],
			vcf
			)
		versions = versions.mix(JVARKIT_VCFGNOMAD.out.versions)
		vcf = JVARKIT_VCFGNOMAD.out.vcf

		if(params.cases!=null && params.controls!=null) {
			CONTRAST1(
				[[id:"cases"],file(params.cases)],
				[[id:"controls"],file(params.controls)],
				vcf
				)
			versions = versions.mix(CONTRAST1.out.versions)
			vcf = CONTRAST1.out.vcf
			}



		BCFTOOLS_INDEX(vcf)
		versions = versions.mix(BCFTOOLS_INDEX.out.versions)
		vcf = BCFTOOLS_INDEX.out.vcf

		if(parseBoolean(params.with_alphamissense)) {
			ALPHAMISSENSE_DOWNLOAD(PREPARE_ONE_REFERENCE.out.dict)
			versions = versions.mix(ALPHAMISSENSE_DOWNLOAD.out.versions)
			
			ALPHAMISSENSE_ANNOTATE(
				ALPHAMISSENSE_DOWNLOAD.out.output,
				vcf
				)
			versions = versions.mix(ALPHAMISSENSE_ANNOTATE.out.versions)
			vcf = ALPHAMISSENSE_ANNOTATE.out.vcf
			}

		

		BCFTOOLS_CONCAT(
			vcf.map{_meta,vcffile,tbi->[vcffile,tbi]}
				.flatMap()
				.collect()
				.map{files->[metadata,files.sort()]}
			)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)




		if(params.cases!=null && params.controls!=null) {
			JVARKIT_VCF_FLATTEN(BCFTOOLS_CONCAT.out.vcf.map{m,f,t->[m,f]})
			versions = versions.mix(JVARKIT_VCF_FLATTEN.out.versions)
			
			CONTRAST2(
				[[id:"cases"],file(params.cases)],
				[[id:"controls"],file(params.controls)],
				JVARKIT_VCF_FLATTEN.out.vcf
				)
			versions = versions.mix(CONTRAST2.out.versions)
			}


		GROUP_BY_GENE(
			[[id:"cases"],(params.cases==null?file(params.cases):([]))],
			[[id:"controls"],(params.controls==null?file(params.controls):([]))],
			BCFTOOLS_CONCAT.out.vcf.map{m,f,t->[m,f]}
			)
		versions = versions.mix(GROUP_BY_GENE.out.versions)

		JVARKIT_VCFSTATS(
			[[id:"sample2pop"],[]],
			BCFTOOLS_CONCAT.out.vcf
			)
		versions = versions.mix(JVARKIT_VCFSTATS.out.versions)
		multiqc = multiqc.mix(JVARKIT_VCFSTATS.out.json.flatMap{flatMapByIndex(it,1)})

		COMPILE_VERSIONS(versions.collect())
		
		MULTIQC(
			[[id:"noconfig"],[]] ,
			multiqc
				.map{meta,f->f}
				.mix(COMPILE_VERSIONS.out.multiqc)
				.collect()
				.map{files->[metadata,files.sort()]}
			)

	}



