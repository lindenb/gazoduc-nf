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
/** https://github.com/FINNGEN/regenie-pipelines */
// https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf/-/blob/d9f89be7c2cccc39336ea282bd54e26146ca645e/workflows/regenie/main.nf

include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { runOnComplete                            } from '../../../modules/utils/functions.nf'
include { parseBoolean                             } from '../../../modules/utils/functions.nf'
include { makeKey                                  } from '../../../modules/utils/functions.nf'
include { isBlank                                  } from '../../../modules/utils/functions.nf'
include { flatMapByIndex                           } from '../../../modules/utils/functions.nf'
include { REGENIE_STEP1                            } from '../../../modules/regenie/step1'
include { DIGEST_SAMPLESHEET                       } from '../step2/sub.nf'
include { MDS_TO_COVARIATES                        } from '../step2/sub.nf'
include { UPDATE_PGEN                              } from '../step2/sub.nf'


workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()
	metadata = [id:"regenie_step1"]
	
	if(params.bim_bed_fam==null) {
		log.error("--bim_bed_fam undefined")
		exit -1
		}
	else if(params.bim_bed_fam.endsWith(".bed")) {
		bgen_ch = Channel.of("${params.bim_bed_fam}")
			.map{F->F.endsWith(".bed")?F.substring(0,F.length()-4):F}
			.map{F->[
				[id:file(F).name],
				file(F+".bed"),
				file(F+".bim"),
				file(F+".fam")
				]
			}
		}
	else
		{
		bgen_ch = Channel.of("${params.bim_bed_fam}")
			.map{F->[
				[id:file(F).name],
				file(F+".bed"),
				file(F+".bim"),
				file(F+".fam")
				]
			}
		}
	



	/** extract pedigrees, cases from vcf input and samplesheet */
	DIGEST_SAMPLESHEET (
		bgen_ch.map{meta,bed,bim,fam->[meta,fam]},
		[[id:"samplesheet"],file(params.samplesheet)]
		)
	versions = versions.mix(DIGEST_SAMPLESHEET.out.versions)


	keep_markers_ch = [[id:"pvar"],(params.step1_markers==null?([]):file("${params.step1_markers}"))]


	/** update pgen with data from pedigee */
	/*
	UPDATE_PGEN(
		DIGEST_SAMPLESHEET.out.plink_ped,
		PLINK2_MERGE_PGEN.out.pgen
		)
	versions = versions.mix(UPDATE_PGEN.out.versions)
	*/

	 /**
	 * COVARIATES, convert from MDS if needed
	 */
	if(params.covariates==null) {
		log.error("--covariates undefined")
		exit -1;
		}
	else if(params.covariates.endsWith(".mds")) {
		MDS_TO_COVARIATES([[id:"covariates"],file(params.covariates)])
		versions = versions.mix(MDS_TO_COVARIATES.out.versions)
		covar_ch = MDS_TO_COVARIATES.out.covariates
	} else {
		covar_ch = [[id:"covariates"],file(params.covariates)]
	}

	
	REGENIE_STEP1( 
		covar_ch /* covariates */,
		DIGEST_SAMPLESHEET.out.plink_ped, 
		keep_markers_ch,
		[[id:"no_exclude_markers"],[]],
		bgen_ch /* bgen file */
		)
	loco_ch = REGENIE_STEP1.out.loco
		.map{meta,locos->[meta,(locos instanceof List ? locos:[locos])]}
		.flatMap{row->flatMapByIndex(row,1)}
	versions = versions.mix(REGENIE_STEP1.out.versions)	
	}

runOnComplete(workflow)
