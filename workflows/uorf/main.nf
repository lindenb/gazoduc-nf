nextflow.enable.dsl=2
include {UORF} from '../../subworkflows/uorf'

workflow {
	genome = Channel.of(file(params.fasta), file(params.fai), file(params.dict)).collect();
	UORF(genome,Channel.fromPath(params.vcf),file(params.bed),file(params.samples))
	}
