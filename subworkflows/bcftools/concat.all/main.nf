include {BCFTOOLS_CONCAT_CONTIGS} from '../concat.contigs'

def ALL_CONTIGS="ALL"


workflow BCFTOOLS_CONCAT_ALL {
take:
	vcfs // tuple [vcf,idx]
	bed //optional bed file
main:
	ch2 = vcfs.
		map{[ALL_CONTIGS,it[0],it[1]]}

	ch4 = BCFTOOLS_CONCAT_CONTIGS(ch2, bed)
emit:	
	output = ch4.output

}

