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

include { runOnComplete                            } from '../../modules/utils/functions.nf'
include { parseBoolean                             } from '../../modules/utils/functions.nf'
include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { PLOT_COVERAGE_01                         } from '../../subworkflows/plotdepth'
include { WALLY                                    } from '../../subworkflows/wally'
include { GTF_INPUT                                } from '../../subworkflows/nf/gtf_input'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { COMPILE_VERSIONS                         } from '../../modules/versions'
include { MULTIQC                                  } from '../../modules/multiqc'
include { BED_STATS                                } from '../../modules/jvarkit/bedstats'

workflow {
 	versions = Channel.empty()
	multiqc = Channel.empty()
	metadata = [id:"plotcoverage"]

	
	if(params.fasta==null) {
		log.warn("--fasta missing")
		exit -1
		}
	if(params.bed==null) {
		log.warn("--bed missing")
		exit -1
	}
	if(params.samplesheet==null) {
		log.warn("--samplesheet missing")
		exit -1
	}

	/***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata.plus(skip_scatter:true),
			Channel.of(params.fasta).map{f->file(f)}.map{f->[[id:f.baseName],f]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	
	GTF_INPUT(
		metadata.plus(
			arg_name:"gtf",
			require_index:true,
			download:true,
			path: params.gtf
		),
		PREPARE_ONE_REFERENCE.out.dict
		)
	versions = versions.mix(GTF_INPUT.out.versions)

	def bed =      [ [id:"userbed"], file(params. bed)  ]
	BED_STATS(bed)
	versions = versions.mix(BED_STATS.out.versions)
	multiqc = multiqc.mix(BED_STATS.out.multiqc)


    READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        )
    versions = versions.mix(READ_SAMPLESHEET.out.versions)

    META_TO_BAMS(
                metadata,
                PREPARE_ONE_REFERENCE.out.fasta,
                PREPARE_ONE_REFERENCE.out.fai,
                READ_SAMPLESHEET.out.samplesheet
                )
    versions = versions.mix(META_TO_BAMS.out.versions)

	PLOT_COVERAGE_01(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
		GTF_INPUT.out.gtf,
		bed,
		META_TO_BAMS.out.bams
		)
	versions = versions.mix(PLOT_COVERAGE_01.out.versions)
	multiqc = multiqc.mix(PLOT_COVERAGE_01.out.multiqc)

	if(parseBoolean(params.with_wally)) {
		WALLY(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			Channel.of(bed).map{meta,bed_file->bed_file}
				.splitCsv(sep:'\t',header:false)
				.map{row->[row[0],(row[1] as int)+1 ,(row[2] as int)]}
				.map{row->{
					def len =  row[2]-row[1];
					def factor = 2.0;
					def len2 = (len*factor+1) as int;
					return [contig:row[0], start: (row[1]-len2>0?row[1]-len2:1), end: row[2]+len2];
					}}
				.combine(META_TO_BAMS.out.bams)
				.map{row,meta,bam,bai->[meta.plus(row),bam,bai]}
			)
		versions = versions.mix(WALLY.out.versions)
		multiqc = multiqc.mix(WALLY.out.multiqc)
		}

	COMPILE_VERSIONS(versions.collect().map{it.sort()})

	MULTIQC(
		[[id:"noconfig"],[]],
		multiqc.map{meta,f->f}
			.mix(COMPILE_VERSIONS.out.multiqc)
			.collect()
			.map{files->[[id:"plotcoverage"],files.sort()]}
		)

}

runOnComplete(workflow);

