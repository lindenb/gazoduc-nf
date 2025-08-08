include {SPLIT_N_VARIANTS  as SPLIT_N_VAR         } from '../../../modules/jvarkit/splitnvariants/main.nf'

List toFlatMap(v) {
	def L1=[];
	def L2 = v[1];
	if(!(L2 instanceof List)) L2=[L2]
	for(f in L2) {
		L1.add([v[0],f])
		}
	return L1;
	}

workflow SPLIT_N_VARIANTS {
take:
    meta
    bed //must be defined, use SCATTER to bed if needed
    vcfs
main:
    versions = Channel.empty()
    

	SPLIT_N_VAR(
		vcfs.combine(bed
			).map{[it[0],it[1],it[2],it[4]]}
		)
	versions= versions.mix(SPLIT_N_VAR.out.versions)

	ch1 = SPLIT_N_VAR.out.vcf.
		flatMap{toFlatMap(it)}

	ch2= SPLIT_N_VAR.out.tbi
		.flatMap{toFlatMap(it)}
	
	vcf = ch1.combine(ch2)
		.filter{it[0].equals(it[2])}
		.filter{(it[1].toRealPath()+".tbi").equals(it[3].toRealPath())}
		.map{[it[0],it[1],it[3]]}
emit:
    versions
    vcf
}
