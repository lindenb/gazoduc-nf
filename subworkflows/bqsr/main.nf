include {GATK4_BASE_RECALIBRATOR} from '../../modules/gatk/base.recalibrator'
include {GATK4_GATHER_BQSR} from '../../modules/gatk/gatherbqsr'
include {GATK4_APPLY_BQSR} from '../../modules/gatk/applybqsr'

workflow BQSR {
    take:
        meta
        fasta
        fai
        dict
        known_vcf // meta,vcf,vcf_idx
        bam
    main:
        versions= Channel.empty()

        beds_ch = CONTIGS_IN_BAM(fasta,fai,dict,bam)
		versions = versions.mix(beds_ch.versions)

		GATK4_BASE_RECALIBRATOR(
            fasta,
            fai,
            dict,
            known_vcf
            CONTIGS_IN_BAM.output.bed.flatMap(T->T[0].map{[T[1],T[2],T[3],it]})
            )
		versions = versions.mix(GATK4_BASE_RECALIBRATOR.out.versions)
	
		GATK4_GATHER_BQSR(GATK4_BASE_RECALIBRATOR.out.table.groupTuple())
		versions = versions.mix(GATK4_GATHER_BQSR.out.versions)

		GATK4_APPLY_BQSR(
            fasta,
            fai,
            dict,
            bam.join(GATK4_GATHER_BQSR.out.table)
            )
		versions = versions.mix(GATK4_APPLY_BQSR.out.versions)
    emit:
        versions
        bam = GATK4_APPLY_BQSR.out.bam
        
}



process CONTIGS_IN_BAM {
tag "${meta.id?:bam.name}"
array 100
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
    tuple val(meta2),path(dict)
	tuple path("BEDS/*.bed"),val(meta),path(bam),path(bai)
output:
	path("bam.contigs.tsv"),emit:output
	path("versions.yml"),emit:versions
script:
    def bqsr_cluster_method="todo"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP BEDS

samtools idxstats --threads ${task.cpus} "${bam}" |\\
    awk -F '\t' '(\$3!=0 && \$1!="*") {printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
    jvarkit  bedcluster.jar --reference "${fasta}" -o BEDS  ${meta.bqsr_cluster_method?:""}
 
cat << EOF > versions.yml
EOF
"""
}
