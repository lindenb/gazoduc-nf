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
nextflow.enable.dsl=2


include { PREPARE_ONE_REFERENCE               } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                        } from '../../subworkflows/samtools/meta2bams1'
include { isBlank                             } from '../../modules/utils/functions.nf'
include { runOnComplete                       } from '../../modules/utils/functions.nf'
include { EXPANSION_HUNTER_DE_NOVO            } from '../../subworkflows/expansion.hunter/denovo'



if( params.help ) {
    exit 0
    }
else
	{
	}



workflow {
	def metadata =[
		id:"xhunterdenovo"
		]
	versions = Channel.empty()
	multiqc = Channel.empty()
	def workflow_meta = [
		id: "pihat"
		]
	if(params.fasta==null) {
		throw new IllegalArgumentException("undefined --fasta");
		}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}

	PREPARE_ONE_REFERENCE(
		metadata.plus([skip_scatter:true]),
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	ch1 = Channel.fromPath(params.samplesheet).splitCsv(header:true, sep:',')

	ch1.filter{it.status.equals("case")}
		.count()
		.filter{it==0}
		.map{throw new IllegalArgumentException("no status=case for ${params.samplesheet}")}
	
	ch1.filter{it.status.equals("control")}
		.count()
		.filter{it==0}
		.map{throw new IllegalArgumentException("no status=control for ${params.samplesheet}")}


	ch1.filter{isBlank(it.status) || !it.status.matches("case|control")}
		.map{throw new IllegalArgumentException("no/invalid status defined for ${it}")}
	
	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		ch1
		)
  	versions = versions.mix(META_TO_BAMS.out.versions)

	EXPANSION_HUNTER_DE_NOVO(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		META_TO_BAMS.out.bams
		)
	
	}

runOnComplete(workflow);
