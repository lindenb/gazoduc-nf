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
include {PREPARE_REFERENCE          } from '../../subworkflows/samtools/prepare.ref'
include {ORA_TO_FASTQ               } from '../../subworkflows/ora/ora2fastq'
include {BAM_TO_FASTQ               } from '../../modules/samtools/bam2fastq'
include {META_TO_PED                } from '../../subworkflows/pedigree/meta2ped'
include {MULTIQC                    } from '../../subworkflows/multiqc'
include {BAM_QC                     } from '../../subworkflows/bamqc'
include {SCATTER_TO_BED             } from '../../subworkflows/gatk/scatterintervals2bed'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


boolean hasKey(def h, def id) {
	return h!=null && h[id]!=null && !(h[id].trim().isEmpty() || h[id].equals("."));
	}

Map cleanupHash(Map h) {
	return h.findAll{k,v->!k.matches("fasta|fai|dict|bam|bai|ora|fastq_1|fastq_2|bed")}
	}

workflow {
  def hash_ref= [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
      ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ hash_ref, file(params.fasta)]
	
	versions = Channel.empty()
	multiqc_ch = Channel.empty()
	ch0 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{assertKeyExistsAndNotEmpty(it,"sample")}
        .map{h->hasKey(h,"id")?h:h.plus(id:h.sample)}
        
     META_TO_PED(hash_ref,ch0)
     versions = versions.mix(META_TO_PED.out.versions)
        
     ch0 = ch0.branch {
        	fastq :  hasKey(it,"fastq_1") && !hasKey(it,"ora") && !hasKey(it,"bam")
        	ora   : !hasKey(it,"fastq_1") &&  hasKey(it,"ora") && !hasKey(it,"bam")
        	bam   : !hasKey(it,"fastq_1") && !hasKey(it,"ora") &&  hasKey(it,"bam") && hasKey(it,"fasta")
        	other : true
        	}
     
     ch0.other.map{throw new IllegalArgumentException("undefined input ${it}.");}
     
     
     ORA_TO_FASTQ(
     	Channel.of(hash_ref),
     	ch0.ora.map{[cleanupHash(it),file(it.ora)]}
     	)
     versions = versions.mix(ORA_TO_FASTQ.out.versions)
     
     BAM_TO_FASTQ(
     	ch0.bam.map{[
     		cleanupHash(it),
     		file(it.bam),
     		(hasKey(it,"bai")?file(it.bai): file(it.bam+(it.bam.endsWith(".bam")?".bai":".crai"))),
     		file(it.fasta),
     		(hasKey(it,"fai")?file(it.fai): file(it.fasta+".fai")),
     		(hasKey(it,"bed")?file(it.bed): [] )
     		]}
     	)
     versions = versions.mix(BAM_TO_FASTQ.out.versions)
     
     ch1 = ch0.fastq
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
        
    
    ch2a = ch1.paired
    	.map{assertKeyExistsAndNotEmpty(it,"fastq_2")}
    	.map{[
		    cleanupHash(it),
		    file(it.fastq_1),
		    file(it.fastq_2)
		    ]}

    ch2b = ch1.single.map{[
        cleanupHash(it),
        file(it.fastq_1),
        []
        ]}

	bed = Channel.of([[id:"nobed"],[]]).first()
	
	PREPARE_REFERENCE(hash_ref,fasta)
	versions = versions.mix(PREPARE_REFERENCE.out.versions)
	fai = PREPARE_REFERENCE.out.fai
	dict = PREPARE_REFERENCE.out.dict
	
	if(params.bwa_index_directory!=null) {
		BWADir = [[id:"bwaindex"],file(params.bwa_index_directory)];
		}
	else
		{
		BWA_INDEX(fasta)
		versions = versions.mix(BWA_INDEX.out.versions)
		BWADir = BWA_INDEX.out.bwa_index
		}


	vcf_for_bqsr= Channel.empty()

	if(params.known_sites!=null) {
		vcf_for_bqsr = Channel.of([hash_ref,[ file(params.known_sites), file(params.known_sites+".tbi")] ])
		}
	else
		{
		vcf_for_bqsr = Channel.of([[id:"novcf"],[] ])
		}


	MAP_BWA(
		hash_ref.plus(
			with_bqsr: (params.known_sites==null || params.with_bqsr==false?false:true),
			with_cram : params.with_cram,
			with_markdup: params.with_markdup,
			markdup_method : params.markdup_method,
			with_seqkit_split : params.with_seqkit_split
			),
		fasta,
		fai,
		dict,
		BWADir,
		vcf_for_bqsr.first(),
		bed,
		ch2a
			.mix(ch2b)
			.mix(ORA_TO_FASTQ.out.fastqs)
			.mix(BAM_TO_FASTQ.out.fastq.flatMap{[
				[it[0],it[1],it[2]],
				[it[0],it[3],[]],
				[it[0],it[4],[]]
				]})
		)
	versions = versions.mix(MAP_BWA.out.versions)
	multiqc_ch = multiqc_ch.mix(MAP_BWA.out.multiqc)

	if(params.capture==null) {
		SCATTER_TO_BED(hash_ref,fasta,fai,dict)
		versions = versions.mix(SCATTER_TO_BED.out.versions)
		bed4qc = SCATTER_TO_BED.out.bed
	  } else {
		bed4qc = Channel.of([hash_ref,file(params.capture)])
	  }

	
	MAKE_SAMPLESHEET(MAP_BWA.out.crams.ifEmpty(MAP_BWA.out.bams)
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
		hash_ref,
		fasta,
		fai,
		dict,
		bed4qc,
		MAP_BWA.out.crams.ifEmpty(MAP_BWA.out.bams)
		)

	MULTIQC(
		hash_ref.plus("id":"bwa"),
		META_TO_PED.out.sample2collection,
		versions,
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

stub:
"""
touch versions.yml samplesheet.csv
"""
}
