/*

Copyright (c) 2023 Pierre Lindenbaum

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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference();

gazoduc.make("bams","NO_FILE").
        description("File containing the paths to the BAM/CRAMS files. One path per line").
        required().
        existingFile().
        put()

gazoduc.make("bed","NO_FILE"). /* not bedS because graphtyper split input by 50kb */
        description("Call in that bed file.").
        required().
        existingFile().
        put()

gazoduc.make("sample2depth","NO_FILE").
        description("tab delimited SAMPLE / DEPTH . Optional. Prevents depth calculation").
        put()

gazoduc.make("dbsnp","NO_FILE").
        description("dbsnp vcf file").
        put()


gazoduc.make("mapq",-1).
        description("mapping quality").
        setInteger().
        put()


if( params.help ) {
    gazoduc.usage().
        name("graphtyper genotype").
        desc("genotype bams using graphtyper").
        print();
    exit 0
} else {
   gazoduc.validate();
}




include {GRAPHTYPER_GENOTYPE_BAMS_01} from '../../../subworkflows/graphtyper/graphtyper.genotype.bams.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete} from '../../../modules/utils/functions.nf'


workflow {
	c1_ch = GRAPHTYPER_GENOTYPE_BAMS_01(params,params.reference,file(params.bams),file(params.bed), file(params.sample2depth) , file(params.dbsnp))
	//html = VERSION_TO_HTML(params,c1_ch.version)	
	}


runOnComplete(workflow);

