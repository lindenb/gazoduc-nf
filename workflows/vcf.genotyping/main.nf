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
include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { BCFTOOLS_MERGE                           } from '../../modules/bcftools/merge3'
include { BCFTOOLS_CONCAT                          } from '../../modules/bcftools/concat3'
include { BCFTOOLS_GENOTYPE                        } from '../../modules/bcftools/genotype'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams2'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { GATK_GENOTYPE_VCF                        } from '../../modules/gatk/genotypevcf'
include { makeKey                                  } from '../../modules/utils/functions.nf'
include { SPLIT_N_VARIANTS                         } from '../../modules/jvarkit/splitnvariants'
include { flatMapByIndex                           } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { BCFTOOLS_STATS                           } from '../../modules/bcftools/stats'
include { MULTIQC                                  } from '../../modules/multiqc'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { JVARKIT_VCF_SET_DICTIONARY as VCFSETDICT1} from '../../modules/jvarkit/vcfsetdict'
include { JVARKIT_VCF_SET_DICTIONARY as VCFSETDICT2} from '../../modules/jvarkit/vcfsetdict'
include { GRAPHTYPER                               } from '../../modules/graphtyper/genotype'
include { MINIGENOTYPER                            } from '../../modules/jvarkit/minigenotyper'
include { BCFTOOLS_QUERY                           } from '../../modules/bcftools/query'
include { BEDTOOLS_SLOP                            } from '../../modules/bedtools/slop'
include { BEDTOOLS_MERGE                           } from '../../modules/bedtools/merge'
include { BIM_TO_VCF                               } from '../../modules/plink/bim2vcf'

List makeArray(meta, List array0,int n) {
	def L = [];
	def i=0;
	while(i< array0.size()) {
		def L2=[];
		while(i< array0.size() && L2.size()<n) {
			L2.add(array0[i]);
			i++;
			}
		L.add([meta,L2]);
		}
	return L;
	}

workflow {
		versions =Channel.empty()
		multiqc = Channel.empty()

		if(!workflow.stubRun) {
			validateParameters()
			}

		if( params.help ) {
			log.info(paramsHelp())
			exit 0
			} 
		else {
			// Print summary of supplied parameters
			log.info paramsSummaryLog(workflow)
			}

		if(params.fasta==null) {
			log.error("undefined --fasta")
			exit -1
			}
		//def main_fasta_id = file(params.fasta).toRealPath().toString().md5()

		if(params.samplesheet==null) {
			log.error("undefined --samplesheet")
			exit -1
			}

		if(params.vcf==null) {
			log.error("undefined --vcf")
			exit -1
			}

		def metadata = [id:"genotyping"]

		PREPARE_ONE_REFERENCE(
			metadata.plus(skip_scatter:true),
			Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
			)
		versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	/* Read samplesheet */
		READ_SAMPLESHEET(
			metadata.plus([arg_name:"samplesheet"]),
			params.samplesheet
			)
		versions = versions.mix(READ_SAMPLESHEET.out.versions)

		/** extract BAM/bai/fasta/fai/dict from samplesheet */
		META_TO_BAMS(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			READ_SAMPLESHEET.out.samplesheet
			)
		versions = versions.mix(META_TO_BAMS.out.versions)


		//add fasta-id
		bams = META_TO_BAMS.out.bams.map{meta,bam,bai,fasta,fai,dict->[
			meta.plus([fasta_id:fasta.toRealPath().toString().md5()]),
			bam,
			bai,
			fasta,
			fai,
			dict
			]}

		all_references = bams.map{meta,bam,bai,fasta,fai,dict->[
				[id:meta.fasta_id,fasta_id: meta.fasta_id],
				fasta,
				fai,
				dict
				]}
			.groupTuple()
			.map{meta,fastas,fais,dicts->[meta,fastas.sort()[0],fais.sort()[0],dicts.sort()[0]]}
		
		all_references.count().map{c->"number of references: ${c}"}.view()

		bams= bams.map{meta,bam,bai,fasta,fai,dict->[meta,bam,bai]}

		if(params.vcf.endsWith(".bim")) {
			bim_ch = [[id:"bim"],file(params.vcf)]
			BIM_TO_VCF(
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
				bim_ch
				)
			versions = versions.mix(BIM_TO_VCF.out.versions)
			vcfs = BIM_TO_VCF.out.vcf
			}
		else
			{
			VCF_INPUT(metadata.plus([
				path: params.vcf,
				arg_name: "vcf",
				require_index : true,
				required: true,
				unique : false
				]))
			versions = versions.mix(VCF_INPUT.out.versions)

			
			vcfs = VCF_INPUT.out.vcf
				.map{meta,vcf,tbi->[meta.plus(id:vcf.toRealPath().toString().md5()),vcf,tbi]}
			
			}

		/* split variants in small chunks */
		if((params.n_variants as int)>0) {
			SPLIT_N_VARIANTS(
				[[id:"nobed"],[]],
				vcfs
				)
			versions = versions.mix(SPLIT_N_VARIANTS.out.versions)
			
			vcfs = SPLIT_N_VARIANTS.out.vcf.flatMap(row->flatMapByIndex(row,1))
				.combine(SPLIT_N_VARIANTS.out.tbi.flatMap(row->flatMapByIndex(row,1)))
				.filter{meta1,vcf,meta2,tbi->meta1.id==meta2.id && "${vcf.name}.tbi" == tbi.name}
				.map{meta1,vcf,meta2,tbi->[meta1.plus(id:vcf.name.md5()),vcf,tbi]}
			}
		
		/** convert VCF input for each type of reference */
		dispatch_ch = vcfs.combine(all_references)
			.multiMap{meta1,vcf,tbi,meta2,fasta,fai,dict->
				dict: [meta2, dict]
				vcf: [meta1.plus(fasta_id:meta2.fasta_id),vcf,tbi]
				}
		
		VCFSETDICT1(
			dispatch_ch.dict,
			dispatch_ch.vcf
			)
		versions = versions.mix(VCFSETDICT1.out.versions)
		

		/** group BAMS by n_bams bams and by meta.fasta_id*/
		grouped_ch = bams
			.map{meta,bam,bai->[[fasta_id:meta.fasta_id],[bam,bai]]}
			.groupTuple()
			.map{meta,arrays->[meta,arrays.sort{a,b->a[0]<=>b[0]}]}
			.flatMap{meta,array->makeArray(meta, array, (params.n_bams as int))}
			.map{meta,array->[meta.plus([id:makeKey(array)]),array.flatten().sort()]}
			
		if(params.method.equals("gatk")) {
			bamvcf_ch = grouped_ch
				.combine(VCFSETDICT1.out.vcf)
				.filter{meta1,bam_array,meta2,vcf,tbi->meta1.fasta_id==meta2.fasta_id}
				.combine(all_references)
				.filter{meta1,bam_array,meta2,vcf,tbi,meta3,fasta,fai,dict->meta1.fasta_id==meta3.fasta_id}
				.multiMap{meta1,bams,meta2,vcf,tbi,meta3,fasta,fai,dict->
					fasta: [meta3,fasta]
					fai: [meta3,fai]
					dict: [meta3,dict]
					vcf : [meta2,vcf,tbi]
					bam : [meta1.plus(vcf_id:meta2.id),bams]
				}
			
			GATK_GENOTYPE_VCF(
				bamvcf_ch.fasta,
				bamvcf_ch.fai,
				bamvcf_ch.dict,
				[[id:"nodbsnp"],[],[]],//no dbsnp because may use multiple fasta
				bamvcf_ch.vcf,
				bamvcf_ch.bam
				)
			versions = versions.mix(GATK_GENOTYPE_VCF.out.versions)
			vcf = GATK_GENOTYPE_VCF.out.vcf
			}
		else if(params.method.equals("graphtyper")) {

			BCFTOOLS_QUERY(VCFSETDICT1.out.vcf)
			versions = versions.mix(BCFTOOLS_QUERY.out.versions)

			BEDTOOLS_SLOP(
				PREPARE_ONE_REFERENCE.out.fai,
				BCFTOOLS_QUERY.out.output
				)
			versions = versions.mix(BEDTOOLS_SLOP.out.versions)

			BEDTOOLS_MERGE(BEDTOOLS_SLOP.out.bed)
			versions = versions.mix(BEDTOOLS_MERGE.out.versions)

			bamvcf_ch = grouped_ch
				.combine(VCFSETDICT1.out.vcf.join(BEDTOOLS_MERGE.out.bed))
				.filter{meta1,bam_array,meta2,vcf,tbi,bed->meta1.fasta_id==meta2.fasta_id}
				.combine(all_references)
				.filter{meta1,bam_array,meta2,vcf,tbi,bed,meta3,fasta,fai,dict->meta1.fasta_id==meta3.fasta_id}
				.multiMap{meta1,bams,meta2,vcf,tbi,bed,meta3,fasta,fai,dict->
					fasta: [meta3,fasta]
					fai: [meta3,fai]
					bed : [meta2,bed]
					vcf : [meta2,vcf,tbi]
					bam : [meta1.plus(vcf_id:meta2.id),bams]
				}


			GRAPHTYPER(
				bamvcf_ch.fasta,
				bamvcf_ch.fai,
				[[id:"nocovlen"],[]],
				bamvcf_ch.vcf,
				bamvcf_ch.bed,
				bamvcf_ch.bam
				)
			versions = versions.mix(GRAPHTYPER.out.versions)
			vcf = GRAPHTYPER.out.vcf
			}
		else if(params.method.equals("jvarkit")) {
			

			bamvcf_ch = grouped_ch
				.combine(VCFSETDICT1.out.vcf)
				.filter{meta1,bam_array,meta2,vcf,tbi->meta1.fasta_id==meta2.fasta_id}
				.combine(all_references)
				.filter{meta1,bam_array,meta2,vcf,tbi,meta3,fasta,fai,dict->meta1.fasta_id==meta3.fasta_id}
				.multiMap{meta1,bams,meta2,vcf,tbi,meta3,fasta,fai,dict->
					fasta: [meta3,fasta]
					fai: [meta3,fai]
					dict: [meta3,dict]
					vcf : [meta2,vcf]
					bam : [meta1.plus(vcf_id:meta2.id),bams]
				}
			
			MINIGENOTYPER(
				bamvcf_ch.fasta,
				bamvcf_ch.fai,
				bamvcf_ch.dict,
				bamvcf_ch.vcf,
				bamvcf_ch.bam
				)
			versions = versions.mix(MINIGENOTYPER.out.versions)
			vcf = MINIGENOTYPER.out.vcf
			}
		else if(params.method.equals("bcftools")) {
			
			VCF2TABIX(VCFSETDICT1.out.vcf.map{meta,vcf,tbi->[meta,vcf]})
			versions = versions.mix(VCF2TABIX.out.versions)


			bamtabix_ch = grouped_ch.combine(VCF2TABIX.out.tabix)
				.filter{meta1,bam_array,meta2,tabix,tbi->meta1.fasta_id==meta2.fasta_id}
				.combine(all_references)
				.filter{meta1,bam_array,meta2,tabix,tbi,meta3,fasta,fai,dict->meta1.fasta_id==meta3.fasta_id}
				.multiMap{meta1,bams,meta2,tabix,tbi,meta3,fasta,fai,dict->
					fasta: [meta3,fasta]
					fai: [meta3,fai]
					tabix : [meta2,tabix,tbi]
					bam : [meta1.plus(vcf_id:meta2.id),bams]
				}

			BCFTOOLS_GENOTYPE(
				bamtabix_ch.fasta,
				bamtabix_ch.fai,
				bamtabix_ch.tabix,
				bamtabix_ch.bam
				)
			versions = versions.mix(BCFTOOLS_GENOTYPE.out.versions)
			vcf = BCFTOOLS_GENOTYPE.out.vcf
			}
		else
			{
			exit(-1,"unknown method ${params.method}");
			}

		
		VCFSETDICT2(
			PREPARE_ONE_REFERENCE.out.dict,
			vcf
			)
		versions = versions.mix(VCFSETDICT2.out.versions)


		BCFTOOLS_MERGE(
			VCFSETDICT2.out.vcf.map{meta,vcf,tbi->[[id:meta.vcf_id],[vcf,tbi]]}
				.groupTuple()
				.map{meta,files->[meta,files.flatten().sort()]}
			)
		versions = versions.mix(BCFTOOLS_MERGE.out.versions)

		BCFTOOLS_CONCAT(
			BCFTOOLS_MERGE.out.vcf.map{meta,vcf,tbi->[vcf,tbi]}
				.flatMap()
				.collect()
				.map{files->[[id:"genotyping"],files.sort()]}
			)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

		if( params.with_stats==true) {
			BCFTOOLS_STATS(
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				[[id:"nobed"],[]],
				[[id:"gtf"],[]],
				[[id:"samples"],[]],
				BCFTOOLS_CONCAT.out.vcf.map{meta,vcf,tbi->[meta,[vcf,tbi]]}
				)
			versions = versions.mix(BCFTOOLS_STATS.out.versions)
			multiqc = multiqc.mix(BCFTOOLS_STATS.out.stats)
			}
		if( params.with_multiqc == true) {
			MULTIQC(
				[[id:"no_mqc_config"],[]],
				multiqc.map{it[1]}.collect().map{[[id:"genotyping"],it]}
				)
			}
	}



process VCF2TABIX {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.02.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(vcf)
output:
	tuple val(meta),path("*.tsv.gz"),path("*.tsv.gz.tbi"),emit:tabix
	path("versions.yml"),emit:versions
script:
	def prefix = task.prefix?:"${meta.id}"
"""
mkdir -p TMP
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\\n' '${vcf}' |\\
	sort -t '\t' -k1,1 -k2,2n -S ${task.memory.kilo} -T TMP |\\
	bgzip > TMP/jeter.tsv.gz

tabix -f -s 1 -b 2 -e 2 TMP/jeter.tsv.gz

mv TMP/jeter.tsv.gz "${prefix}.tsv.gz"
mv TMP/jeter.tsv.gz.tbi "${prefix}.tsv.gz.tbi"
touch versions.yml
"""

stub:
	def prefix = task.prefix?:"${meta.id}"
"""
touch "${prefix}.tsv.gz" "${prefix}.tsv.gz.tbi" versions.yml
"""
}
