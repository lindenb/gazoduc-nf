/*

Copyright (c) 2025 Pierre Lindenbaum

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

include {MULTIQC  as MQC1            } from '../../modules/multiqc'
include {MULTIQC  as MQC2            } from '../../modules/multiqc'
include {COMPILE_VERSIONS           } from '../../modules/versions'
include {JVARKIT_MULTIQCPOSTPROC    } from '../../modules/jvarkit/multiqcpostproc'

workflow MULTIQC {
take:
	meta
	sample2collection
	versions 
	multiqc_ch
main:
	COMPILE_VERSIONS(versions.collect().map{it.sort()})
    multiqc_ch = multiqc_ch.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MQC1(multiqc_ch.map{it[1]}.collect().map{[meta,it]})
	
    JVARKIT_MULTIQCPOSTPROC(sample2collection.filter{meta,sample2pop->(sample2pop?true:false)}, MQC1.out.datadir)
    
    MQC2(JVARKIT_MULTIQCPOSTPROC.out.json.map{m,json->[m.plus("id":m.id+".sample2pop"),json]})
    

}