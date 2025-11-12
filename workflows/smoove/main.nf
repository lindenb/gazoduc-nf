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

include {SMOOVE                        } from '../../subworkflows/smoove' 
include {SMOOVE_ANNOTATE               } from '../../modules/smoove/annotate' 
include {JVARKIT_VCF_SET_DICTIONARY    } from '../../modules/jvarkit/vcfsetdict' 
include {runOnComplete                 } from '../../modules/utils/functions.nf'
include {isBlank                       } from '../../modules/utils/functions.nf'
include {parseBoolean                  } from '../../modules/utils/functions.nf'
include {META_TO_BAMS                  } from '../../subworkflows/samtools/meta2bams1'
include {PREPARE_ONE_REFERENCE         } from '../../subworkflows/samtools/prepare.one.ref'
include {GFF3_INPUT                    } from '../../subworkflows/nf/gff3_input'
include {READ_SAMPLESHEET              } from '../../subworkflows/nf/read_samplesheet'
include {COMPILE_VERSIONS              } from '../../modules/versions'


workflow {
	if(params.fasta==null) {
		throw new IllegalArgumentException("undefined --fasta");
		}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}

	 versions = Channel.empty()
    def metadata = [
            id: "smoove"
            ]
     
	def exclude    = [[id:"no_exclude"], [] ]


	
	PREPARE_ONE_REFERENCE(
		metadata, 
		Channel.fromPath(params.fasta).map{[[id:it.simpleName],it]}
		)
	fasta= PREPARE_ONE_REFERENCE.out.fasta.first()
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	/* no fastq samplesheet */
    samplesheet0_ch = READ_SAMPLESHEET(
        [arg_name:"samplesheet"],
        params.samplesheet
        ).samplesheet
	versions = versions.mix(READ_SAMPLESHEET.out.versions)


	

	if(params.exclude_bed==null) {
		GET_EXCLUDE(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			PREPARE_ONE_REFERENCE.out.complement_bed
			)
		versions  = versions.mix(GET_EXCLUDE.out.versions)
		exclude = GET_EXCLUDE.out.bed
		}
	else {
		exclude    = [[id:"user_xclude"], file(params.exclude_bed) ]
		}
	
	META_TO_BAMS(
        metadata ,
        PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        samplesheet0_ch
        )
	bams =  META_TO_BAMS.out.bams
		.map{meta,bam,bai->[meta.plus("status":(isBlank(meta.status)?"case":meta.status)),bam,bai]}
	versions  = versions.mix(META_TO_BAMS.out.versions)

	SMOOVE(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
        PREPARE_ONE_REFERENCE.out.fai,
        PREPARE_ONE_REFERENCE.out.dict,
		exclude,
		bams
		)
	versions  = versions.mix(SMOOVE.out.versions)

	vcf = SMOOVE.out.vcf
	if(parseBoolean(params.with_annotation)) {
		if(params.gff3==null) {
			log.err("undefined --gff3 for annotation")
			exit -1
			}
			
		GFF3_INPUT(
				[
				arg_name:"gff3",
				path: "${params.gff3}",
				require_index :true,
				download : true
				],
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict
			)
		versions = versions.mix(GFF3_INPUT.out.versions)
		
		SMOOVE_ANNOTATE(
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			GFF3_INPUT.out.gff3,
			SMOOVE.out.vcf
			)
		versions  = versions.mix(SMOOVE_ANNOTATE.out.versions)
		vcf = SMOOVE_ANNOTATE.out.vcf
		}
	
	
	JVARKIT_VCF_SET_DICTIONARY(
		PREPARE_ONE_REFERENCE.out.dict,
		vcf
		)
	versions  = versions.mix(JVARKIT_VCF_SET_DICTIONARY.out.versions)
	vcf = JVARKIT_VCF_SET_DICTIONARY.out.vcf
	

	COMPILE_VERSIONS(versions.collect().map{it.sort()})
	}



process GET_EXCLUDE {
tag "${fasta.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
	tuple val(meta4),path(complement_bed)
output:
	tuple val(meta1),path("*.bed"),emit:bed /* exclude2.bed otherwise collision with bed from scatter_bed */
    path("versions.yml"),emit:versions
script:
 def prefix = task.ext.prefix?:"${meta1.id}.exclude"
"""
hostname 1>&2

mkdir -p TMP

#cf https://github.com/brentp/smoove

if ${meta1.ucsc_name=="hg19"}
then
	echo 'https://github.com/hall-lab/speedseq/blob/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed' > TMP/jeter.url
fi

if ${meta1.ucsc_name=="hg38"}
then
	echo 'https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed'  > TMP/jeter.url
fi

if test -f TMP/jeter.url
then
	curl -L `cat TMP/jeter.url` |\\
	cut -f 1,2,3 |\\
	jvarkit bedrenamechr -R "${fasta}" --column 1 --convert SKIP > TMP/exclude2.bed 

fi

awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {printf("%s\t0\t%s\\n",\$1,\$2);}' '${fai}' >> TMP/exclude2.bed
cat ${complement_bed} >> TMP/exclude2.bed

sort -T TMP -t '\t' -k1,1 -k2,2n TMP/exclude2.bed > ${prefix}.bed

test -s ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
    url: \$(cat TMP/jeter.url || true)
END_VERSIONS
"""


stub:
def prefix = task.ext.prefix?:"${meta1.id}.exclude"

"""
touch versions.yml ${prefix}.bed
"""
}


runOnComplete(workflow)
