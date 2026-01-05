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
include {UNMAPPED                                 } from '../../subworkflows/unmapped'
include {runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { SAMTOOLS_STATS                          } from '../../modules/samtools/stats'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include {PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include {META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams2'
include { READ_SAMPLESHEET                        } from '../../subworkflows/nf/read_samplesheet'


workflow {
	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
			throw new IllegalArgumentException("undefined --samplesheet");
			}



	def versions = Channel.empty()
	def multiqc  = Channel.empty()
	def metadata = [id:"contaminations"]
	def acn_file = [metadata,file("${moduleDir}/contaminants.txt")];
	if(params.acn_file!=null) {
		acn_file = [[id:file(params.acn_file).name],file(params.acn_file)]
		}	

  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata,
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/* no fastq samplesheet */
	READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        )
	versions = versions.mix(READ_SAMPLESHEET.out.versions)

	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)

	


	UNMAPPED(
		metadata,
		acn_file,
		META_TO_BAMS.out.bams
		)
	versions = versions.mix(UNMAPPED.out.versions)
	multiqc = multiqc.mix(UNMAPPED.out.multiqc)

 	COMPILE_VERSIONS(versions.collect())
	multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

    MULTIQC(
		[[id:"no_mqc_config"],[]],
		multiqc.collect().map{[[id:"contaminations"],it]}
		)
	}


runOnComplete(workflow);
