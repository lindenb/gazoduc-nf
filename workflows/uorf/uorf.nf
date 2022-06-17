nextflow.enable.dsl=2
include {UORF} from '../../subworkflows/uorf/uorf.nf'
workflow {
	UORF(params,params.reference,params.vcf,params.bed)
	}
