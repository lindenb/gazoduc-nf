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

include {dumpParams;runOnComplete   } from '../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty   } from '../../modules/utils/functions.nf'
include {assertKeyExistsAndNotEmpty } from '../../modules/utils/functions.nf'
include {FASTQC                     } from '../../modules/fastqc'
include {MULTIQC                    } from '../../subworkflows/multiqc'
include {META_TO_PED                } from '../../subworkflows/pedigree/meta2ped'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}




workflow {
	if(params.samplesheet==null) {
		log.error("undefined --samplesheet");
		exit -1;
		}

	versions = Channel.empty()
	multiqc_ch = Channel.empty()
	ch0 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true,sep:",")
        
    META_TO_PED([id:"fastqc"],ch0)
    versions = versions.mix(META_TO_PED.out.versions)    
        
        
    ch1 = ch0.map{assertKeyExistsAndNotEmpty(it,"sample")}
        .map{
            if(it.fastq_1==null && it.R1!=null) return it.plus(fastq_1:it.R1);
            return it;
            }
        .map{
            if(it.fastq_2==null && it.R2!=null) return it.plus(fastq_2:it.R2);
            return it;
            }  
        .map{assertKeyExistsAndNotEmpty(it,"fastq_1")}
        .branch{
                paired: it.fastq_2!=null &&  !it.fastq_2.isEmpty() &&  !it.fastq_2.equals(".")
                single: true
                }
    
    
    
    
    ch2a = ch1.paired.map{assertKeyExistsAndNotEmpty(it,"fastq_2")}.map{[
        [id:it.sample],
        [file(it.fastq_1),file(it.fastq_2)]
        ]}

    ch2b = ch1.single.map{[
        [id:it.sample],
        [file(it.fastq_1)]
        ]}
	
	

    FASTQC( ch2a.mix(ch2b) )
    versions = versions.mix(FASTQC.out.versions)
    multiqc_ch = multiqc_ch.mix(FASTQC.out.zip)
    multiqc_ch = multiqc_ch.mix(FASTQC.out.html)

    MULTIQC(
		["id":"multiqc"],
		META_TO_PED.out.sample2collection,
		versions,
        [[id:"no_mqc_config"],[]],
		multiqc_ch
		)
    }


runOnComplete(workflow);

