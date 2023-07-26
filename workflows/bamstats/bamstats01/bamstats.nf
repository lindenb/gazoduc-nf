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


def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults()


gazoduc.build("bams","NO_FILE").
        desc("list of BAM/CRAM files").
        required().
        existingFile().
        put()

gazoduc.build("sample2pop","NO_FILE").
        desc("optional path to TSV : sample(TAB)pop to plot boxplot").
        put()


include {isBlank;runOnComplete} from '../../../modules/utils/functions.nf'
include {SAMTOOLS_STATS_01} from '../../../subworkflows/samtools/samtools.stats.01.nf'


if( params.help ) {
    gazoduc.usage().
	name("SAMTOOLS stats").
	desc("samtools stats").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch1 = SAMTOOLS_STATS_01([:],file(params.bams), file(params.sample2pop) )
	//PUBLISH(ch1.zip)
	}

process PUBLISH {
tag "${zip.name}"
input:
	path(zip)
when:
	!isBlank(params.getOrDefault("publishDir",""))
script:
"""
mkdir -p "${params.publishDir}"
cp -v "${zip}" "${params.publishDir}/"
"""
}


runOnComplete(workflow);
