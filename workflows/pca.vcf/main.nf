/**

Thank you Floriane Simonet for the Help

*/


include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include {moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {ACP_VCF_STEP as ACP_STEP01} from '../../subworkflows/pca.vcf/step.pca.vcf.01.nf' addParams([step_id:"pca_step1",step_name:"Step 1"])
include {MULTIQC} from '../../subworkflows/multiqc/multiqc.nf'


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)






runOnComplete(workflow)

workflow {
	genome_ch = Channel.value([file(params.fasta),file(params.fasta+".fai"),file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict")])
	ACP_VCF(genome_ch, file(params.vcf), file(params.sample2collection), file(params.blacklisted_bed),file(params.exclude_samples) )
	}

workflow ACP_VCF {
	take:
		genome_ch
		vcf
		sample2collection
		blacklisted_bed
		samples
	main:
		step1_ch = ACP_STEP01(genome_ch, vcf, sample2collection, blacklisted_bed,samples)

		to_multiqc = step1_ch.multiqc


                mqc_ch = MULTIQC(to_multiqc)

	}


