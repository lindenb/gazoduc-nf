/*

Copyright (c) 2026 Pierre Lindenbaum

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
include { MULTIQC                     } from '../../../modules/multiqc'
include { BCFTOOLS_STATS              } from '../../../modules/bcftools/stats'
include { BCFTOOLS_CALL               } from '../../../subworkflows/bcftools/call'
include { COMPILE_VERSIONS            } from '../../../modules/versions/main.nf'
include { runOnComplete; dumpParams   } from '../../../modules/utils/functions.nf'
include { JVARKIT_BAM_RENAME_CONTIGS  } from '../../../modules/jvarkit/bamrenamechr'
include { BEDTOOLS_MAKEWINDOWS        } from '../../../modules/bedtools/makewindows'
include { BED_CLUSTER                 } from '../../../modules/jvarkit/bedcluster'
include { PREPARE_ONE_REFERENCE       } from '../../../subworkflows/samtools/prepare.one.ref'
include { PREPARE_USER_BED            } from '../../../subworkflows/bedtools/prepare.user.bed'
include { META_TO_BAMS                } from '../../../subworkflows/samtools/meta2bams2'
include { META_TO_PED                 } from '../../../subworkflows/pedigree/meta2ped'
include {isBlank                      } from '../../../modules/utils/functions'



if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}




workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()

	if(params.fasta==null) {
		throw new IllegalArgumentException("undefined --fasta");
		}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}

	

	def metadata=[
		id: "bcftools"
		]



	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta)
			.map{[[id:it.baseName],it]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	def gtf  =    [ metadata, file(params.gtf) ]
	def pedigree = [ metadata, []]

	
	if(params.bed==null) {
		bed = Channel.of([ [id:"nobed"], [] ])
		}
	else {
		PREPARE_USER_BED(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			PREPARE_ONE_REFERENCE.out.scatter_bed,
			Channel.of([[id:"capture"],file(params.bed)])
			)
		versions = versions.mix(PREPARE_USER_BED.out.versions)
		multiqc = multiqc.mix(PREPARE_USER_BED.out.multiqc)
		bed = PREPARE_USER_BED.out.bed
		}
	
	if(params.samplesheet.endsWith(".json")) {
		ch0 = Channel.fromPath(params.samplesheet)
        	.splitJSon()
		}
	else
		{
		ch0 = Channel.fromPath(params.samplesheet)
        	.splitCsv(header:true,sep:',')
		}

	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		ch0
		)
	versions = versions.mix(META_TO_BAMS.out.versions)


	META_TO_PED(metadata, META_TO_BAMS.out.bams.map{it[0]})
	versions = versions.mix(META_TO_PED.out.versions)


	bams_ch = META_TO_BAMS.out.bams
        .map{
			if(!isBlank(it[0].batch)) return it;
			def L=[ it.plus(batch:it.id) ];
			L.addAll(it.subsList(1,it.size()));
			return L;
			}
		.branch {
			ok_ref: java.nio.file.Files.isSameFile(it[3],file(params.fasta))
			bad_ref: true
			}
		

	JVARKIT_BAM_RENAME_CONTIGS(
		PREPARE_ONE_REFERENCE.out.dict,
		bed,
		bams_ch.bad_ref
		)

  versions =versions.mix(JVARKIT_BAM_RENAME_CONTIGS.out.versions)
		
  BEDTOOLS_MAKEWINDOWS(bed)
  versions =versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)
  bed = BEDTOOLS_MAKEWINDOWS.out.bed
   
  BED_CLUSTER(
	PREPARE_ONE_REFERENCE.out.dict,
	bed
	)
  versions =versions.mix(BED_CLUSTER.out.versions)

  bed = BED_CLUSTER.out.bed.flatMap{
		def L=[];
		for(f in it[1]) {
			L.add([[id:f.name],f]);
			}
		return L;
		}

	BCFTOOLS_CALL(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
        META_TO_PED.out.pedigree,
		bed,
		bams_ch.ok_ref
			.map{meta,bam,bai,_fasta,_fai,_dict->[meta,bam,bai]}
			.mix(JVARKIT_BAM_RENAME_CONTIGS.out.bam)
		)
	versions = versions.mix(BCFTOOLS_CALL.out.versions)
	



	BCFTOOLS_STATS(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		bed,
		gtf,
		[[id:"no_samples"],[]],
		BCFTOOLS_CALL.out.vcf.map{[it[0],[it[1],it[2]]]}
		)
	versions = versions.mix(BCFTOOLS_STATS.out.versions)
	multiqc = multiqc.mix(BCFTOOLS_STATS.out.stats.map{it[1]})
	

	COMPILE_VERSIONS(versions.collect())
	multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)
	
	MULTIQC(
		[[id:"no_mqc_config"],[]],
		multiqc.collect().map{[[id:"bcftools"],it]}
		)
	}

runOnComplete(workflow)


