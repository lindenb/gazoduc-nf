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

/*
	sn_ch = pedigree.
		.splitText(elem:1)
		.map{[it[0],it.trim().split("[ \t]+")]}
		.filter{!it[1].startsWith("#)}
		.map{[it[0].plus([id:it[1][1], father:it[1][2], mother:it[1][3]])]} //meta
		.filter{!(it[0].father.equals("0") || it[0].mother.equals("0"))}
	

	trio3ch = sn2.sn_ch
		.combine(bams) //meta1,C.meta,C.bam,C.bai
		.filter{it[0].id.equals(it[1].id)}
		.map{[it[0].id,it[1],it[2]]} // meta,C.bam,C.bai
		.combine(bams) //  meta,C.bam,C.bai, F.meta, F.bam, F.bai
		.filter{it[0].father.equals(it[3].id)}
		.map{[it[0],it[1],it[2],it[4],it[5]]} // meta,C.bam,C.bai,F.bam,F.bai
		.combine(bams) //  meta,C.bam,C.bai, F.meta, F.bam, F.bai, M.meta, M.bam,M.bai
		.filter{it[0].mother.equals(it[6].id)}
		.map{[it[0],it[1],it[2],it[3],it[4],it[5],it[7],it[8]]} // meta,C.bam,C.bai,F.bam,F.bai
	
	// find child without two parent and try to find one of the two parent
	trio2_ch = sn_ch.map{[it[0],null]} // add null as a placeholder to flag child that cannot be join-ed
		.join(trio3ch)  ,remainder:true)
		.filter{it[1]==null} // child has no father AND mother
		.map{[it[0]]}
		.combine(bams) //meta1,C.meta,C.bam,C.bai
		.filter{it[0].id.equals(it[1].id)}
		.map{[it[0].id,it[1],it[2]]} // meta,C.bam,C.bai
		.combine(bams) //  meta,C.bam,C.bai, F.meta, F.bam, F.bai
		.filter{it[0].father.equals(it[3].id || it[0].mother.equals(it[3].id)}
		.map{[it[0].plus(parent:it[3].id),it[1],it[2],it[4],it[5]]} // meta,C.bam,C.bai,P.bam,P.bai
	
	used_ch = trio3ch.flatMap{[
			[id:it[0].id] , 
			[id:it[0].father]
			[id:it[0].mother],
			]}.mix(flatMap{[
				[id:it[0].id],
				[id:it[0].parent]
				]}.map{[it]}
	
	// find other unmapped samples
	no_trio_ch = bams
		.map{[it[0].id,it[0],it[1],it[2]]}
		.join(used_ch,remainder:true)
		.filter{it[1]==null}
		.map{[it[1],it[2],it[3]]} //meta, bam,bai

	gvcfs= Channel.empty()

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
	

	DEEP_TRIO3(
		fasta,
		fai,
		dict,
		trio2ch.map{it.plus(bed[1)}.plus(bed[1)}
		)
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
	
	
emit:
	versions
	vcf = GLNEXUS_GENOTYPE.out.vcf
*/
vcf= Channel.empty()
versions
}
