/*

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include {BED_CLUSTER          } from '../../../modules/jvarkit/bedcluster'
include {DEEP_TRIO3           } from '../../../modules/deepvariant/deep_trio3'
include {DEEP_VARIANT_CALL    } from '../../../modules/deepvariant/call'
include {GLNEXUS_GENOTYPE     } from '../../../modules/glnexus/genotype'
include {BCFTOOLS_CONCAT       } from '../../../modules/bcftools/concat'


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
	no_trio_ch = Channel.empty()

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

		no_trio_ch = Channel.empty()//TODO a ce jour pas de situation ou les echantillons ne sont pas utilises
		}
	else	{
		no_trio_ch = bams
		}

	DEEP_TRIO3(
		fasta,
		fai,
		dict,
		trio3_ch.combine(bed)
			.map{[it[0],it[2],it[3],it[5],it[6],it[8],it[9],it[11]]} //meta,c.bam,c.bai, f.bam,f.bai, m.bam,m.bai, bed
		)
	
	versions= versions.mix(DEEP_TRIO3.out.versions)


	DEEP_VARIANT_CALL(
		fasta,
		fai,
		dict,
		no_trio_ch.combine(bed)
			.map{[it[0],it[1],it[2],it[4]]}.view() //meta,bam,bai,bed
		)
	versions= versions.mix(DEEP_VARIANT_CALL.out.versions)

	to_genotype = Channel.empty()
	//TODO add DEEP_VARIANT_CALL, TRIO2


	to_genotype = to_genotype.mix(
		DEEP_TRIO3.out.gvcf
			.map{[
				it[3].toRealPath(),/*bed*/
				it[1].plus(it[2]) /* vcf and index */
				]}
		)
	/** group by bed */
	to_genotype  = to_genotype
		.groupTuple()
		.map{[it[0],it[1].flatten()]}
		

	ch1 = to_genotype.multiMap{
			bed : [[itd:file(it[0]).baseName], it[0]]
			gvcf :[[itd:file(it[0]).baseName], it[1]]
		}
	
	
	GLNEXUS_GENOTYPE(
		ch1.bed,
		[[id:"no_config"],[]],
		ch1.gvcf
		)
	versions= versions.mix(GLNEXUS_GENOTYPE.out.versions)
	

	BCFTOOLS_CONCAT(
		GLNEXUS_GENOTYPE.out.vcf
			.map{[it[1],it[2]]} //vcf, tbi
			.collect()
			.map{[[id:"deepvariant"],it]},
		[[id:"no_bed"],[]]
		)
	versions= versions.mix(BCFTOOLS_CONCAT.out.versions)

emit:
	vcf = BCFTOOLS_CONCAT.out.vcf
	versions
}
