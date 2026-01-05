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
include { runOnComplete                 } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE         } from '../../../subworkflows/samtools/prepare.one.ref'
include { MULTIQC                       } from '../../../modules/multiqc'
include { CONCORDANCE_BAMS              } from '../../../subworkflows/concordance/bams1'
include { COMPILE_VERSIONS              } from '../../../modules/versions/main.nf'

workflow {
	if(params.samplesheet==null) {
		log.warn("undefined --samplesheet");
		java.lang.System.exit(-1);
		}
	if(params.fasta==null) {
		log.warn("undefined --fasta");
		java.lang.System.exit(-1);
		}
	if(params.vcf==null) {
		log.warn("undefined --vcf");
		java.lang.System.exit(-1);
		}
    
   versions = Channel.empty()
   multiqc = Channel.empty()
   metadata = [
		id:"concordance"
		]

    PREPARE_ONE_REFERENCE(
   		[id:"x",skip_scatter:true],
   		Channel.of(params.fasta).map{f->file(f)}.map{f->[[id:f.baseName],file(f)]}.first()
   		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

    bams = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true,sep:',')
        .map{row->[
                [
                    id:row.sample,
                    sample:row.sample
                ],
                file(row.bam),
                file(row.bam+".crai")
                ]}
    vcf = Channel.of([[id:file(params.vcf).baseName],file(params.vcf),file(params.vcf+".tbi")])
    
    CONCORDANCE_BAMS(
        metadata,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
        vcf,
        bams
        )
    versions = versions.mix(CONCORDANCE_BAMS.out.versions)

    COMPILE_VERSIONS(versions.collect())
    multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

    MULTIQC(
                [[id:"no_mqc_config"],[]],
                multiqc.collect().map{[[id:metadata.id],it]}
                )

    }

runOnComplete(workflow)
