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

gazoduc.build("bams", "NO_FILE").
	desc(gazoduc.Gazoduc.DESC_BAM_OR_CRAM_LIST).
	existingFile().
	required().
	put()


gazoduc.build("reference_in", "NO_FILE").
	desc("Reference for input BAMS. "+ gazoduc.Gazoduc.DESC_INDEXED_FASTA).
	existingFile().
	required().
	indexedFasta().
	put()

gazoduc.build("reference_out", "NO_FILE").
	desc("Reference for output BAMS. "+ gazoduc.Gazoduc.DESC_INDEXED_FASTA).
	existingFile().
	required().
	indexedFasta().
	put()



include {REMAP_BWA_01} from '../../../subworkflows/mapping/remap.bwa.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {runOnComplete} from '../../../modules/utils/functions.nf'



if( params.help ) {
    gazoduc.usage().
	name("Remap bwa").
	desc("Takes a mapped set of BAM or CRAM files and remap them on a new reference").
	print();
    exit 0
} else {
   gazoduc.validate();
}


workflow {
	remap_ch = REMAP_BWA_01(params,params.reference_in,params.reference_out,params.bams)
	//html = VERSION_TO_HTML(params,remap_ch.version)	
	}

runOnComplete(workflow);

