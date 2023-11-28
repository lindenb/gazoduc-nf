/**

Thank you Floriane Simonet for the Help

*/

include {moduleLoad;runOnComplete;dumpParams} from '../../modules/utils/functions.nf'
include {ACP_VCF_STEP as ACP_STEP01} from '../../subworkflows/pca.vcf/step.pca.vcf.01.nf'
include {MULTIQC} from '../../subworkflows/multiqc/multiqc.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


runOnComplete(workflow)

workflow {
	ACP_VCF([:], params.genomeId, file(params.vcf), file(params.sample2collection))
	}

workflow ACP_VCF {
	take:
		meta
		genomeId
		vcf
		sample2collection
	main:
		step1_ch = ACP_STEP01([:],genomeId, vcf, sample2collection)

		to_multiqc = step1_ch.multiqc


                mqc_ch = MULTIQC(to_multiqc)

	}

