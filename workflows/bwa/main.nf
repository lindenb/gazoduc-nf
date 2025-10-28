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
include {MAP_BWA                    } from '../../subworkflows/bwa/map.fastqs'
include {BWA_INDEX                  } from '../../modules/bwa/index'
include {runOnComplete              } from '../../modules/utils/functions.nf'
include {PREPARE_ONE_REFERENCE      } from '../../subworkflows/samtools/prepare.one.ref'
include {META_TO_PED                } from '../../subworkflows/pedigree/meta2ped'
include {MULTIQC                    } from '../../subworkflows/multiqc'
include {BAM_QC                     } from '../../subworkflows/bamqc'
include {IF_EMPTY                   } from '../../subworkflows/nf/if_empty'
include {SAMPLESHEET_TO_FASTQ       } from '../../subworkflows/samplesheet2fastq'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


boolean hasKey(def h, def id) {
	return h!=null && h[id]!=null && !(h[id].trim().isEmpty() || h[id].equals("."));
	}

workflow {

	if(params.samplesheet==null) {
		System.err.println("undefined --samplesheet");
		System.exit(-1);
		}
	if(params.fasta==null) {
		System.err.println("undefined --fasta");
		System.exit(-1);
		}

  	def workflow_medadata = [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
      ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ workflow_medadata, file(params.fasta)]
	
	versions = Channel.empty()
	multiqc_ch = Channel.empty()


	/**
	 *
	 * LOAD SAMPLESHEET
	 *
	 */

	ch0 =  Channel.fromPath(params.samplesheet);

	if(params.samplesheet.endsWith(".json")) {
		ch0 = ch0.splitJson()
		}
	else
		{
		ch0 =  ch0.splitCsv(header:true,sep:',')
		}

	ch0 = ch0
		.map{assertKeyExistsAndNotEmpty(it,"sample")}
        .map{h->hasKey(h,"id")?h:h.plus(id:h.sample)}
    
	SAMPLESHEET_TO_FASTQ(
		workflow_medadata.plus([
			bam2fastq_method : params.bam2fastq_method
			]),
		ch0
		)
	versions = versions.mix(SAMPLESHEET_TO_FASTQ.out.versions)

    META_TO_PED(workflow_medadata,ch0)
    versions = versions.mix(META_TO_PED.out.versions)
    
	bed = Channel.of([[id:"nobed"],[]]).first()
	
	PREPARE_ONE_REFERENCE(workflow_medadata,Channel.of(fasta))
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	
	
	if(params.bwa_index_directory!=null) {
		BWADir = [[id:"bwaindex"],file(params.bwa_index_directory)];
		}
	else
		{
		BWA_INDEX(PREPARE_ONE_REFERENCE.out.fasta)
		versions = versions.mix(BWA_INDEX.out.versions)
		BWADir = BWA_INDEX.out.bwa_index
		}


	vcf_for_bqsr= Channel.empty()

	if(params.known_sites!=null) {
		vcf_for_bqsr = Channel.of([
			workflow_medadata,
			[ file(params.known_sites), file(params.known_sites+".tbi")]
			])
		}
	else
		{
		vcf_for_bqsr = Channel.of([[id:"novcf"],[] ])
		}


	MAP_BWA(
		workflow_medadata.plus(
			with_bqsr: (params.known_sites==null || params.with_bqsr==false?false:true),
			with_cram : params.with_cram,
			with_markdup: params.with_markdup,
			markdup_method : params.markdup_method,
			with_seqkit_split : params.with_seqkit_split
			),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		BWADir,
		vcf_for_bqsr.first(),
		bed,
		SAMPLESHEET_TO_FASTQ.out.paired_end
		)

	versions = versions.mix(MAP_BWA.out.versions)
	multiqc_ch = multiqc_ch.mix(MAP_BWA.out.multiqc)

	if(params.capture==null) {
		bed4qc = bed = PREPARE_ONE_REFERENCE.out.scatter_bed.first()
	  } else {
		bed4qc = Channel.of([workflow_medadata,file(params.capture)])
	  }

	bams_out = MAP_BWA.out.bams
	if(params.with_cram) {
		bams_out = MAP_BWA.out.crams
		}
	
	MAKE_SAMPLESHEET(
		bams_out
		.map{meta,bam,bai->[
			meta.id,
			"${params.outdir}/BAMS/${meta.id}/${params.prefix?:""}${bam.name}",
			"${params.outdir}/BAMS/${meta.id}/${params.prefix?:""}${bai.name}",
			meta.sex?:"",
			meta.father?:"",
			meta.mother?:"",
			meta.status?:"",
			meta.collection?:"",
			"${params.fasta}"
			]}
		.map{it.join(",")}
		.collect()
		)
	versions = versions.mix(MAKE_SAMPLESHEET.out.versions)


	
	BAM_QC(
		workflow_medadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		bed4qc,
		IF_EMPTY(MAP_BWA.out.crams,MAP_BWA.out.bams)
		)

	versions = versions.mix(BAM_QC.out.versions)
	multiqc_ch = multiqc_ch.mix(BAM_QC.out.multiqc)

	MULTIQC(
		workflow_medadata.plus("id":"bwa"),
		META_TO_PED.out.sample2collection,
		versions,
		[[id:"no_mqc_config"],[]],
		multiqc_ch
		)
		
    }


runOnComplete(workflow);




process MAKE_SAMPLESHEET {
tag "N=${L.size()}"
input:
	val(L)
output:
	path("samplesheet.csv"),emit:samplesheet
	path("versions.yml"),emit:versions
script:
"""

echo 'sample,bam,bai,sex,father,mother,status,collection,fasta' > jeter.csv
cat << EOF >> jeter.csv
${L.join("\n")}
EOF

mv jeter.csv samplesheet.csv

touch versions.yml
"""
}
