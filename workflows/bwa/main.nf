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
nextflow.enable.dsl=2


include {dumpParams                 } from '../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty   } from '../../modules/utils/functions.nf'
include {assertKeyExistsAndNotEmpty } from '../../modules/utils/functions.nf'
include {FASTQC                     } from '../../modules/fastqc'
include {MULTIQC                    } from '../../modules/multiqc'
include {COMPILE_VERSIONS           } from '../../modules/versions'
include {MAP_BWA                    } from '../../subworkflows/bwa/map.fastqs'
include {BWA_INDEX                  } from '../../modules/bwa/index'
include {runOnComplete              } from '../../modules/utils/functions.nf'
include {PREPARE_REFERENCE          } from '../../subworkflows/samtools/prepare.ref'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}





workflow {
  def hash_ref= [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
      ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ hash_ref, file(params.fasta)]
	def fai   = [ hash_ref, file(params.fai)]
	def dict  = [ hash_ref, file(params.dict)]

	versions = Channel.empty()
	multiqc_ch = Channel.empty()
	ch1 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{assertKeyExistsAndNotEmpty(it,"sample")}
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
        [id:it.sample, single_end:false, paired_end:true],
        [file(it.fastq_1),file(it.fastq_2)]
        ]}

    ch2b = ch1.single.map{[
        [id:it.sample,single_end:true, paired_end:false],
        [file(it.fastq_1),file([])]
        ]}

	bed = Channel.of([[id:"nobed"],[]])
	
	if(params.bwa_index_directory!=null) {
		BWADir = [[id:"bwaindex"],file(params.bwa_index_directory)];
		}
	else
		{
		BWA_INDEX(fasta)
		versions = versions.mix(BWA_INDEX.out.versions)
		BWADir = BWA_INDEX.out.bwa_index
		}


	PREPARE_REFERENCE(fasta)
	versions = versions.mix(PREPARE_REFERENCE.out.versions)
	fai = PREPARE_REFERENCE.out.fai
	dict = PREPARE_REFERENCE.out.dict

	MAP_BWA(
		hash_ref,
		fasta,
		fai,
		dict,
		BWADir,
		[[id:"bqsr"],[]],
		bed,
		ch2a.mix(ch2b)
		)
	versions = versions.mix(MAP_BWA.out.versions)

	COMPILE_VERSIONS(versions.collect().map{it.sort()})
    multiqc_ch = multiqc_ch.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MULTIQC(multiqc_ch.map{it[1]}.collect().map{[[id:"mapbwa"],it]})
    }


runOnComplete(workflow);

