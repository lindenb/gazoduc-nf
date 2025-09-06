/*

Copyright (c) 2024 Pierre Lindenbaum

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

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include {runOnComplete;} from '../../../modules/utils/functions.nf'


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)



include {SOMALIER_BAMS    } from  '../../../subworkflows/somalier/bams'
include {MULTIQC          } from '../../../subworkflows/multiqc/multiqc.nf'




workflow {
	def hash_ref= [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName
		]
	def fasta = [ hash_ref, file(params.fasta)]
	def fai   = [ hash_ref, file(params.fai)]
	def dict  = [ hash_ref, file(params.dict)]
	def sites    = (params.user_sites==null ? [ hash_ref, []]:[hash_ref,file(params.user_sites)])
	def pedigree = (params.pedigree==null   ? [ hash_ref, []]:[hash_ref,file(params.pedigree  )])

	SOMALIER_BAMS(
		[:],
		Channel.of(fasta),
		fai,
		dict,
		Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.map{
				if(!it.containsKey("sample")) throw new IllegalArgumentException("sample missing");
				if(!it.containsKey("bam")) throw new IllegalArgumentException("bam missing");
				if(!it.containsKey("bai")) throw new IllegalArgumentException("bai missing");
				def key=  [id:it.sample];
				if(it.containsKey("fasta") && !it.fasta.trim().isEmpty()) {
					def fasta2 = it.fasta
					def fai2 = (it.containsKey("fai") && !it.fasta.trim().isEmpty()? it.fai : fasta2 +".fai")
					return [key,file(it.bam),file(it.bai),file(fasta2),file(fai2)]
					}
				return [key,file(it.bam),file(it.bai)]
				},
		pedigree,
		sites
		)
	//multiqc =  MULTIQC(somalier_ch.output.mix(somalier_ch.qc).flatten().collect())
	}


//runOnComplete(workflow)
