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
nextflow.enable.dsl=2

//include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'




// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   //log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
//validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)



include {SOMALIER_BAMS              } from  '../../../subworkflows/somalier/bams'
include {MULTIQC                    } from '../../../subworkflows/multiqc'
include {PREPARE_REFERENCE          } from '../../../subworkflows/samtools/prepare.ref'
include {assertKeyExistsAndNotEmpty } from '../../../modules/utils/functions.nf'
include {runOnComplete              } from '../../../modules/utils/functions.nf'
include {META_TO_PED                } from '../../../subworkflows/pedigree/meta2ped'
include {SOMALIER_FIND_SITES        } from '../../../modules/somalier/find.sites'

boolean hasKey(def h, def id) {
	return h!=null && h[id]!=null && !(h[id].trim().isEmpty() || h[id].equals("."));
	}


workflow {
	def versions = Channel.empty()
	def multiqc = Channel.empty()
	def hash_ref= [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
		ucsc_name : "${params.ucsc_name?:"undefined"}"
		]
	
	
	def fasta = [ hash_ref, file(params.fasta)]
	
	PREPARE_REFERENCE(
		hash_ref.plus(skip_scatter:true), 
		Channel.of(fasta)
		)
	fasta= PREPARE_REFERENCE.out.fasta.first()
    versions = versions.mix(PREPARE_REFERENCE.out.versions)

	
	ch0 = Channel.fromPath(params.samplesheet)
		.splitCsv(header:true,sep:',')
		.map{assertKeyExistsAndNotEmpty(it,"sample")}
		.map{assertKeyExistsAndNotEmpty(it,"bam")}
		.map{
			if(hasKey(it,"bai")) return it;
			return it;
			}
		.map{it.bai?it: (it.bam.endsWith(".bam") ? it.plus(["bai":it.bam+".bai"]):  it.plus(["bai":it.bam+".crai"]))}
		.map{it.id?it:it.plus(id:it.sample)}
		
	
	META_TO_PED(hash_ref, ch0)
    versions = versions.mix(META_TO_PED.out.versions)
		
	if(params.pedigree==null) {
		pedigree = META_TO_PED.out.pedigree_gatk
		}
	else
		{
		pedigree = Channel.of([hash_ref,file(params.pedigree)] )
		}
		
	if(params.population_vcf!=null) {
		SOMALIER_FIND_SITES([hash_ref,file(params.population_vcf),file(population_vcf+".tbi")] )
		versions = versions.mix(SOMALIER_FIND_SITES.out.versions)
		sites = SOMALIER_FIND_SITES.out.vcf
		}	
	else if(params.user_sites!=null) {
		sites = Channel.of([hash_ref,file(params.user_sites)])
		}
	else
		{
		sites = Channel.of([ hash_ref, []] )
		}



	SOMALIER_BAMS(
		hash_ref,
		fasta,
		PREPARE_REFERENCE.out.fai.first(),
		PREPARE_REFERENCE.out.dict.first(),
		ch0.map{
			def key=  [id:it.id];
			if(hasKey(it,"fasta")) {
				def fasta2 = it.fasta
				def fai2 = (hasKey(it,"fai")? it.fai : fasta2 +".fai")
				return [key,file(it.bam),file(it.bai),file(fasta2),file(fai2)]
				}
			return [key,file(it.bam),file(it.bai)]
			},
		pedigree,
		sites
		)
	multiqc = multiqc.mix(SOMALIER_BAMS.out.multiqc)
	
	MULTIQC(
		hash_ref.plus("id":"somalier"),
		META_TO_PED.out.sample2collection,
		versions,
		multiqc
		)
	}


runOnComplete(workflow)
