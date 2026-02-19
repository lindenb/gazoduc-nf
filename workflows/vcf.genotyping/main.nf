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
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { GATK_GENOTYPE_VCF                        } from '../../modules/gatk/genotypevcf'
include { makeKey                                  } from '../../modules/utils/functions.nf'
include { SPLIT_N_VARIANTS                         } from '../../modules/jvarkit/splitnvariants'
include { flatMapByIndex                           } from '../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { BCFTOOLS_STATS                           } from '../../modules/bcftools/stats'
include  {MULTIQC                                  } from '../../modules/multiqc'

List makeArray(array0,int n) {
	def L = [];
	def i=0;
	while(i< array0.size()) {
		def L2=[];
		while(i< array0.size() && L2.size()<n) {
			L2.add(array0[i]);
			i++;
			}
		L.add(L2);
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
			metadata,
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
			READ_SAMPLESHEET.out.samplesheet
			)
		versions = versions.mix(META_TO_BAMS.out.versions)

		SPLIT_N_VARIANTS(
			[[id:"nobed"],[]],
			[[id:"vcf"],file(params.vcf),file(params.vcf+".tbi")]
			)
		versions = versions.mix(SPLIT_N_VARIANTS.out.versions)


	    vcfs = SPLIT_N_VARIANTS.out.vcf.flatMap(row->flatMapByIndex(row,1))
            .combine(SPLIT_N_VARIANTS.out.tbi.flatMap(row->flatMapByIndex(row,1)))
            .filter{meta1,vcf,meta2,tbi->meta1.id==meta2.id && "${vcf.name}.tbi" == tbi.name}
            .map{meta1,vcf,meta2,tbi->[meta1.plus(id:vcf.name.md5()),vcf,tbi]}



		/** group BAMS by n_bams bams */
		grouped_ch = META_TO_BAMS.out.bams
			.map{meta,bam,bai->[bam,bai]}
			.toSortedList{a,b->a[0]<=>b[0]}
			.flatMap(array->makeArray(array, (params.n_bams as int)))
			.map{array->[[id:makeKey(array)],array.flatten().sort()]}
			.view()


		if(params.dbsnp==null) {
			dbsnp = 	[[id:"nodbsnp"],[],[]]
			}
		else
			{
			dbsnp = 	[[id:"dbsnp"],file(params.dbsnp),file(params.dbsnp+".tbi")]
			}

		if(params.method.equals("gatk")) {
			bamvcf_ch = grouped_ch.combine(vcfs)
				.multiMap{meta1,bams,meta2,vcf,tbi->
					vcf : [meta2,vcf,tbi]
					bam : [meta1.plus(vcf_id:meta2.id),bams]
				}


			GATK_GENOTYPE_VCF(
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
				dbsnp,//no dbsnp
				bamvcf_ch.vcf,
				bamvcf_ch.bam
				)
			versions = versions.mix(GATK_GENOTYPE_VCF.out.versions)
			vcf = GATK_GENOTYPE_VCF.out.vcf
			}
		else if(params.method.equals("bcftools")) {
			VCF2TABIX(vcfs.map{meta,vcf,tbi->[meta,vcf]})
			versions = versions.mix(VCF2TABIX.out.versions)

			bamtabix_ch = grouped_ch.combine(VCF2TABIX.out.tabix)
				.multiMap{meta1,bams,meta2,tabix,tbi->
					tabix : [meta2,tabix,tbi]
					bam : [meta1.plus(vcf_id:meta2.id),bams]
				}

			BCFTOOLS_GENOTYPE(
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
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
		BCFTOOLS_MERGE(
			vcf.map{meta,vcf,tbi->[[id:meta.vcf_id],[vcf,tbi]]}
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

		MULTIQC(
			[[id:"no_mqc_config"],[]],
			multiqc.map{it[1]}.collect().map{[[id:"genotyping"],it]}
			)
	}



process VCF2TABIX {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../../conda/bioinfo.02.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(vcf)
output:
	tuple val(meta),path("*.tsv.gz"),path("*.tsv.gz.tbi"),emit:tabix
	path("versions.yml"),emit:versions
script:
	def prefix = task.meta.prefix?:"${meta.id}"
"""
mkdir -p TMP
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\\n' '${vcf}' |\\
	sort -t '\t' -k1,1 -k2,2n -S ${task.memory.kilo} -T TMP |\\
	bgzip > TMP/jeter.tsv.gz

tabix -f -s 1 -b 2 -e 2 TMP/jeter.tsv.gz

mv TMP/jeter.tsv "${prefix}.tsv.gz"
mv TMP/jeter.tsv.tbi "${prefix}.tsv.gz.tbi"
touch versions.yml
"""
stub:
	def prefix = task.meta.prefix?:"${meta.id}"
"""
touch "${prefix}.tsv.gz" "${prefix}.tsv.gz.tbi" versions.yml
"""
}