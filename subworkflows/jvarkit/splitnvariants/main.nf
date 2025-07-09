include {SCATTER_TO_BED                           } from '../../gatk/scatterintervals2bed/main.nf'
include {SPLIT_N_VARIANTS                         } from '../../../modules/jvarkit/splitnvariants/main.nf'
include {BED_CLUSTER                              } from '../../../modules/jvarkit/bedcluster/main.nf'

List toFlatMap(v) {
	def L1=[];
	def L2 = v[1];
	if(!(L2 instanceof List)) L2=[L2]
	for(f in L2) {
		L1.add([v[0],f])
		}
	return L1;
	}

workflow SPLIT_VCF {
take:
    meta
    fasta
    fai
    dict
    bed
    vcfs
main:
    versions = Channel.empty()
    
	if(!bed[1]) {
		SCATTER_TO_BED(meta,fasta,fai,dict)
		bed = SCATTER_TO_BED.out.bed
		versions= versions.mix(SCATTER_TO_BED.out.versions)

		BED_CLUSTER(fasta,fai,dict,bed)
		versions= versions.mix(BED_CLUSTER.out.versions)
		bed = BED_CLUSTER.out.bed.flatMap{toFlatMap(it)}
    	} 

	SPLIT_N_VARIANTS(
		vcfs.combine(bed).map{[it[0],it[1],it[2],it[4]]}
		)
	versions= versions.mix(SPLIT_N_VARIANTS.out.versions)

	ch1 = SPLIT_N_VARIANTS.out.vcf.
		flatMap{toFlatMap(it)}

	ch2= SPLIT_N_VARIANTS.out.tbi
		.flatMap{toFlatMap(it)}
	
	vcf = ch1.combine(ch2)
		.filter{it[0].equals(it[2])}
		.filter{(it[1].toRealPath()+".tbi").equals(it[3].toRealPath())}
		.map{[it[0],it[1],it[3]]}
emit:
    versions
    vcf
}