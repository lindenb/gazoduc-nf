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

include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER1} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER2} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER3} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER4} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER5} from "./sub.nf"
include {DIVIDE_AND_CONQUER as DIVIDE_AND_CONQUER6} from "./sub.nf"
include { JVARKIT_BAM_RENAME_CONTIGS               } from '../../../modules/jvarkit/bamrenamechr'
include { PREPARE_ONE_REFERENCE                    } from '../../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                             } from '../../../subworkflows/samtools/meta2bams2'
include { READ_SAMPLESHEET                         } from '../../../subworkflows/nf/read_samplesheet'
include { BEDTOOLS_MAKEWINDOWS                     } from '../../../modules/bedtools/makewindows'
include { PREPARE_USER_BED                         } from '../../../subworkflows/bedtools/prepare.user.bed'
include { MOSDEPTH                                 } from '../../../modules/mosdepth'
include { runOnComplete                            } from '../../../modules/utils/functions.nf'
include { BED_CLUSTER                              } from '../../../modules/jvarkit/bedcluster'
include { BEDTOOLS_SLOP                            } from '../../../modules/bedtools/slop'
include { READ_LENGTH                              } from '../../../modules/graphtyper/read.length'
include { isSameFai                                } from '../../../modules/utils/functions.nf'

workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()
	def metadata = [
		id: "graphtyper"
		]


	if(!workflow.stubRun) {
		validateParameters()
		}

	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

	/* Read samplesheet */
	READ_SAMPLESHEET(
		metadata.plus([arg_name:"samplesheet"]),
		params.samplesheet
		)
	versions = versions.mix(READ_SAMPLESHEET.out.versions)

	/** extract BAM/bai/fasta/fai/dict from samplesheet */
	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)


  	bams = META_TO_BAMS.out.bams
			.combine(PREPARE_ONE_REFERENCE.out.fai)
			.map {meta1,bam,bai,fasta,fai1,dict,meta2,fai2->
					[
					meta1.plus(ok_ref:isSameFai(fai1,fai2)),
					bam,
					bai,
					fasta,
					fai1
					]
				}
	
	if(params.bed!=null) {
		PREPARE_USER_BED(
            metadata,
            PREPARE_ONE_REFERENCE.out.fasta,
            PREPARE_ONE_REFERENCE.out.fai,
            PREPARE_ONE_REFERENCE.out.dict,
            PREPARE_ONE_REFERENCE.out.scatter_bed,
            Channel.of([[id:file(params.bed).baseName],file(params.bed)])
            )
        versions = versions.mix(PREPARE_USER_BED.out.versions)
        multiqc = multiqc.mix(PREPARE_USER_BED.out.multiqc)
        bed = PREPARE_USER_BED.out.bed.first()
		}
	else
		{
		bed = PREPARE_ONE_REFERENCE.out.scatter_bed
		}

	
	
	BEDTOOLS_MAKEWINDOWS(bed)
	versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

	BED_CLUSTER(
			PREPARE_ONE_REFERENCE.out.dict,
			BEDTOOLS_MAKEWINDOWS.out.bed
			)
	versions = versions.mix(BED_CLUSTER.out.versions)
	beds_ch = BED_CLUSTER.out.bed
			.map{_meta,beds->beds}
			.map{beds->beds instanceof List?beds:[beds]}
			.flatMap()
			.map{bed->[[id:bed.baseName],bed]}


	BEDTOOLS_SLOP(
		PREPARE_ONE_REFERENCE.out.fai,
		bed
		)
	versions = versions.mix(BEDTOOLS_SLOP.out.versions)

	JVARKIT_BAM_RENAME_CONTIGS(
		PREPARE_ONE_REFERENCE.out.dict,
		BEDTOOLS_SLOP.out.bed,
		bams
			.filter{meta,bam,bai,fasta,fai,dict->meta.ok_ref==false}
			.map{meta,bam,bai->[meta.plus(depth:-1),bam,bai]}
		)
	versions =versions.mix(JVARKIT_BAM_RENAME_CONTIGS.out.versions)

	bams = bams.filter{meta,bam,bai,fasta,fai,dict->meta.ok_ref==true}
			.mix(JVARKIT_BAM_RENAME_CONTIGS.out.bam)
x	

	bams2 = bams.ok_ref.branch{meta,bam,bai->
			has_depth : meta.depth!=null && !meta.depth.isEmpty() && meta.depth!="." && (meta.depth as int)>0
			no_depth : true
			}
	
	
	MOSDEPTH(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		bams2.no_depth
			.combine(BEDTOOLS_SLOP.out.bed)
			.map{meta1,bam,bai,meta2,bed->[meta1,bam,bai,bed]}
		)
	versions = versions.mix(MOSDEPTH.out.versions)

	/*

	ch1 = MOSDEPTH.out.summary_txt.splitText().
                map{[it[1],it[2],it[3],it[0].trim()]}.
		mix(bams2.has_depth.map{[it.sample,it.bam,it.bai,it.depth]})
		/*
	readlen_ch = READ_LENGTH(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		ch1
		)
	versions = versions.mix(READ_LENGTH.out.versions)

	READ_LENGTH.out..splitText().map{[it[1],it[0].trim()]}.join(ch1)

	samplesheet_ch = DIGEST_DEPTH(ch2.map{it.join("\t")}.collect())


	
        intervals_ch  = intervals_ch.output.
                splitCsv(header:false,sep:'\t').
                map{[it[0],""+((it[1] as int)+1),it[2]]}


        merge_ch = Channel.empty()

        level1 = DIVIDE_AND_CONQUER1(1, samplesheet_ch.output , intervals_ch )
        merge_ch = merge_ch.mix(level1.ok)

        level2 = DIVIDE_AND_CONQUER2(2, samplesheet_ch.output , level1.failed)
        merge_ch = merge_ch.mix(level2.ok)

        level3 = DIVIDE_AND_CONQUER3(3, samplesheet_ch.output , level2.failed)
        merge_ch = merge_ch.mix(level3.ok)

        level4 = DIVIDE_AND_CONQUER4(4, samplesheet_ch.output , level3.failed)
        merge_ch = merge_ch.mix(level4.ok)

        level5 = DIVIDE_AND_CONQUER5(5, samplesheet_ch.output , level4.failed)
        merge_ch = merge_ch.mix(level5.ok)

        level6 = DIVIDE_AND_CONQUER6(6, samplesheet_ch.output , level5.failed)
        merge_ch = merge_ch.mix(level6.ok)

	SAVE_FAILED(level6.failed.map{it.join("\t")}.collect())
	
        MERGE(merge_ch.map{T->{
			String c = T[0];
			if(!c.matches("(chr)?[0-9XY]+")) {
				c="others";
				}
			return [c,T[1]];
			}}.groupTuple())

        level6.failed.view{"${it} cannot be called"}
	*/
	}


runOnComplete(workflow);


process __MOSDEPTH {
tag "${sample}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/mosdepth.yml"
input:
	path(genome)
	path(bed)
	tuple val(sample),path(bam),path(bai)
output:
	tuple path("${sample}.cov.txt"),val(sample),path(bam),path(bai), emit:output
script:
	def fasta = genome.find{it.name.endsWith("a")}
	def mapq = params.mapq
"""
mkdir -p TMP

# bed for autosomes
awk -F '\t' '(\$1 ~/^(chr)?[0-9]+\$/)' '${bed}' > TMP/jeter.bed
test -s TMP/jeter.bed

mosdepth  \\
 	-t ${task.cpus} \\
	--by TMP/jeter.bed \\
	--no-per-base \\
	--fasta "${fasta}" \\
	--mapq ${mapq} \\
	TMP/output \\
	${bam}

awk -F '\t' '(\$1=="total_region") {print \$4}' TMP/output.mosdepth.summary.txt > "${sample}.cov.txt"
test "${sample}.cov.txt"
"""
}




process DIGEST_DEPTH {
label "process_single"
input:
	val(L)
output:
	path("samplesheet.tsv"),emit:output
script:
"""
cat << EOF > jeter.tsv
${L.join("\n")}
EOF

sort -T . -t '\t' -k1,1 jeter.tsv > samplesheet.tsv
test -s samplesheet.tsv
cut -f1 samplesheet.tsv | sort | uniq -d > dups.txt
test ! -s dups.txt
"""
}


process RDF {
input:
	val(L)
output:
	path("digest.ttl"),emit:output
script:
"""

cat << EOF | awk -F '\t' '{printf("samples:%s rdf:type foaf:Person;\\n\tu:Library libraries:\\"%s\".\\n",\$1,\$1); }' > digest.ttl
${L.join("\n")}
EOF

"""
}


process MERGE {
tag "${contig} N=${L.size()}"
label "process_medium"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
        tuple val(contig),val(L)
output:
        path("${contig}.merged.bcf")
        path("${contig}.merged.bcf.csi")
script:
	def args = "--remove-duplicates"
"""
mkdir -p TMP
set -x
cat << EOF >  TMP/jeter.list
${L.join("\n")}
EOF

SQRT=`awk 'END{X=NR;z=sqrt(X); print (z==int(z)?z:int(z)+1);}' "TMP/jeter.list"`
split -a 9 --additional-suffix=.list --lines=\${SQRT} TMP/jeter.list TMP/chunck.


find TMP/ -type f -name "chunck*.list" | while read F
do
                bcftools concat --threads ${task.cpus} ${args} -a -O b --file-list "\${F}" -o "\${F}.bcf" 
		bcftools index --threads ${task.cpus} -f "\${F}.bcf"
                echo "\${F}.bcf" >> TMP/jeter2.list
done

bcftools concat ${args} --threads ${task.cpus} -a -O b9 --file-list TMP/jeter2.list -o "TMP/jeter.bcf" 
bcftools index --threads ${task.cpus} -f "TMP/jeter.bcf"

mv TMP/jeter.bcf ${contig}.merged.bcf
mv TMP/jeter.bcf.csi ${contig}.merged.bcf.csi
"""
}

process SAVE_FAILED {
executor "local"
input:
	val(L)
output:
	path("failed.txt"),emit:output
script:
"""
cat << EOF | sort > failed.txt
${L.join("\n")}
EOF
"""
}
