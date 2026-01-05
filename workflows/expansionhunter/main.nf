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
include { EXPANSION_HUNTER                    } from '../../subworkflows/expansion.hunter'
include { XHUNTER_DOWNLOAD_CATALOG            } from '../../modules/expansion.hunter/download.catalog'
include { COMPILE_VERSIONS                    } from '../../modules/versions'

if( params.help ) {
    exit 0
    }




workflow {
	def metadata =[
		id:"xhunter"
		]
	versions = Channel.empty()
	multiqc = Channel.empty()
	
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

	if(params.catalog ==null) {
		XHUNTER_DOWNLOAD_CATALOG(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai
			)
		versions = versions.mix(XHUNTER_DOWNLOAD_CATALOG.out.versions)
		catalog = XHUNTER_DOWNLOAD_CATALOG.out.catalog
	} else {
		catalog = Channel.of([[id:"catalog"],file(params.catalog)]).first()
	}

	if(params.samplesheet.endsWith(".json")) {
		ch1 = Channel.fromPath(params.samplesheet).splitJson()
		}
	else
		{
		ch1 = Channel.fromPath(params.samplesheet).splitCsv(header:true, sep:',')
		}
	
	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		ch1
		)
  	versions = versions.mix(META_TO_BAMS.out.versions)

	META_TO_BAMS.out.bams
		.filter{meta,_bam,_bai->isBlank(meta.sex) || !meta.sex.matches("(male|female)")}
		.map{throw new IllegalArgumentException("sex undefined/unknown for ${it}")}

	

	EXPANSION_HUNTER(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		catalog,
		META_TO_BAMS.out.bams
		)
	versions = versions.mix(EXPANSION_HUNTER.out.versions)
	COMPILE_VERSIONS(versions.collect().map{it.sort()})
	}

runOnComplete(workflow);