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




//include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {dumpParams;isBlank;runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT} from '../../../subworkflows/bcftools/bcftools.concat.02.nf'
include {GATK3_HC_MINIKIT} from './step.gatk3.nf'
include {GATK4_HC_MINIKIT} from './step.gatk4.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	ch1 = HC_MINIKIT(params.genomeId, file(params.bams), file(params.beds) )
	}

runOnComplete(workflow);


workflow HC_MINIKIT {
	take:
		genomeId
		bams
		beds
	main:
		version_ch = Channel.empty()
		if(beds.name.equals("NO_FILE")) {
			throw new IllegalArgumentException("--beds undefined");
			}
		else	{
			each_bed  = Channel.fromPath(beds).splitText().map{it.trim()}
			}
		if(params.gatkjar.contains("local.jar")) {
			call_ch = GATK4_HC_MINIKIT(genomeId,bams,each_bed)
			}
		else if(params.gatkjar.contains("GenomeAnalysisTK"))
			{
			call_ch = GATK3_HC_MINIKIT(genomeId,bams,each_bed)
			}
		else
			{
			throw new IllegalArgumentException("bah gatk jar ${params.gatkjar}");
			}
		


		ch2 = BCFTOOLS_CONCAT(["method":"all"], call_ch.output.map{T->[vcf:""+T[0],vcf_index:""+T[1]]},file("NO_FILE"))
		version_ch = version_ch.mix(ch2.version)

        version_ch = MERGE_VERSION("hc-gatk-minikit", version_ch.collect())

/*
emit:
        vcfs = ch2.vcfs
        version= version_ch
*/
	}


