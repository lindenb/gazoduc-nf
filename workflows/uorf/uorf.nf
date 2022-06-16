include {UORF} from ' ../../subworkflows/uorf/uorf.nf'
workflow {
	UORF(params,params.reference,params.bed)
	}
