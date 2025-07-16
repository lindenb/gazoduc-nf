include {BED_CLUSTER} from '../../../modules/jvarkit/bedcluster'


workflow DEEPVARIANT_TRIOS {
take:
	meta
	fasta
	fai
	dict
	bed
	pedigree
	bams //[ meta, bam,bai]
	
main:
	versions = Channel.empty()
	

	BED_CLUSTER(fasta,fai,dict,bed)
	versions = versions.mix(BED_CLUSTER.out.versions)

	bed = BED_CLUSTER.out.bed.flatMap{
		def L=[];
		for(f in it[1]) {
			L.add([[id:f.name],f]);
			}
		return L;
		}

	
	trio3_ch=Channel.empty()
	without_family_ch = bams

	if(pedigree[1]) {
		trios_ch = Channel.of(pedigree[1])
			.splitText()
			.map{it.trim().split("[ \t]+")}
			.filter{!(it[2].equals("0") && it[3].equals("0"))} //may be one parent or two parants
			.map{[id:it[1],father:it[2],mother:it[3]]}
			.map{[it]}
			.unique()



		trio3_ch = trios_ch
			.combine(bams) //meta1,C.meta,C.bam,C.bai
			.filter{it[0].id.equals(it[1].id)}
			.combine(bams) //meta1,C.meta,C.bam,C.bai, P.meta,P.bam,P.bai, 
			.filter{it[0].father.equals(it[4].id)}
			.combine(bams) //meta1,C.meta,C.bam,C.bai, P.meta,P.bam,P.bai,   M.meta,M.bam,M.bai, 
			.filter{it[0].mother.equals(it[7].id)}


		// TODO situation ou il n'y a qu'un seul parent
		
		remains1 = trio3_ch.flatMap{[
			[it[0].id,"used"],
			[it[0].father,"used"],
			[it[0].mother,"used"]
			]}.join(bams.map{[it[0].id,it[0],it[1],it[2]]},failOnDuplicate:true, remainder:true )
			.filter{it[1]==null} // ? a verifier
			.view{"TODODOOOOOTO a ce jour pas de situation ou les echantillons ne sont pas utilises"}

		without_family_ch = Channel.empty()//TODO a ce jour pas de situation ou les echantillons ne sont pas utilises
	}


		DEEP_TRIO3(
		fasta,
		fai,
		dict,
		trio3_ch.map{it.plus(bed[1)}.plus(bed[1)}
		)


/*
	
	DEEP_TRIO2(
		fasta,
		fai,
		dict,
		trio2ch.map{it.plus(bed[1)}
		)
	versions = versions.mix(DEEP_TRIO2.out.versions)
	
	gvcfs = gvcfs.mix(DEEP_TRIO2.out.gvcfs
		.map{[
			it[5].toRealPath(),//bed
			[
			it[1],it[2], //C.gvcf
			it[3],it[4]  //P.gvcf
			]
			]})
	


	versions = versions.mix(DEEP_TRIO3.out.versions)
	
	gvcfs = gvcfs.mix(DEEP_TRIO3.out.gvcfs
		.map{[
			it[7].toRealPath(),//bed
			[
			it[1],it[2], //C.gvcf
			it[3],it[4], //F.gvcf
			it[5],it[6]  //M.gvcf
			]
			]})
	
	
	DEEP_VARIANT_CALL(
		fasta,
		fai,
		dict,
		no_trio_ch.map{it.plus(bed[1)}
		)
	versions = versions.mix(DEEP_TRIO3.out.versions)
	
	gvcfs = gvcfs.mix(DEEP_VARIANT_CALL.out.gvcfs
		.map{[
			it[7].toRealPath(),//bed
			[
			it[1],it[2] //C.gvcf
			]})
	
	
	ch1 = gvcfs.groupTuple()
		.map{[it[0],it[1].flatten()} // [bed, gvcfs]
		.map{[[id:it[0]],it[0],it[1]} // [meta,bed, gvcfs]
		.multiMap{ 
			bed: [it[0],it[1]]
			gvcfs :  [it[0],it[2]]
		}
	
	GLNEXUS_GENOTYPE(
		ch1.bed,
		[[id:"no_config"],[]],
		ch1.gvcfs
		)
	versions = versions.mix(GLNEXUS_GENOTYPE.out.versions)
	
*/
	
emit:
	versions
}
