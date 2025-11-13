/**
ORIGINAL snakemake workflow by Raphael Blanchet PhD.


*/

process PB_GERMLINE {
  tag "${meta.id} ${R1.name} ${R2.name}"
  label 'process_short'
  label 'parabricks'

  // afterScript "rm -rf TMP"
  input:
	path(fasta)
	path(fai)
	path(known_vcf)
	path(known_vcf_tbi)
	tuple val(meta),path(R1),path(R2)
  output:
	tuple val(meta),
		path("${meta.id}.cram"),
		path("${meta.id}.cram.crai"),emit:bams
	tuple val(meta),
		path("${meta.id}.g.vcf.gz"),
		path("${meta.id}.g.vcf.gz.tbi"),
		emit:gvcf , optional:true
        tuple val(meta),
		path("${meta.id}.duplicate.metrics.txt"), emit:bqsr, optional:true
	path("versions.yml"),emit:versions
  script:
    def sample = meta.id
    def lib = meta.LIB?:sample
    def pl = meta.PL?:"ILLUMINA"
    def with_bqsr = task.ext.with_bqsr==true
    def bwa="-M"
    def cpu_per_gpu = task.ext.cpu_per_gpu
    def bqsr_args = with_bqsr? "--knownSites \"${known_vcf.name}\"   --out-recal-file \"${sample}.bqsr.report.txt\"   " : ""
    def fixmate_args = task.ext.with_fixmate ? "--fix-mate" :""
    def low_memory = task.ext.low_memory==true || task.attempt>1 ? "--low-memory" : ""
 """
	mkdir -p TMP

           pbrun germline \\
                --num-gpus ${task.ext.gpus} \\
                --gpusort \\
                --bwa-cpu-thread-pool ${cpu_per_gpu} \\
                --num-htvc-threads ${cpu_per_gpu} \\
                --num-cpu-threads-per-stage ${cpu_per_gpu} \\
                --memory-limit ${task.memory.giga} \\
                --ref ${fasta} \\
                --in-fq "${R1}" "${R2}" \\
                --out-bam "${sample}.cram" \\
                --out-variants "${sample}.g.vcf.gz" \\
                --read-group-sm "${sample}" \\
                --read-group-lb "${lib}" \\
                --read-group-pl "${pl}" \\
                --gpusort \\
		--gvcf \\
		--out-duplicate-metrics "${sample}.duplicate.metrics.txt" \\
		${bqsr_args} \\
		${fixmate_args} \\
                --tmp-dir TMP \\
                --logfile ${sample}.log \\
		${low_memory}


cat <<-END_VERSIONS > versions.yml
"${task.process}":
        pbrun: \$(grep Parabr -m1 -A1 "${sample}.log" | grep -o 'Version [^ ]*' )
END_VERSIONS
"""


    stub:
    def sample = meta.id
    """
    touch ${sample}.cram
    touch ${sample}.cram.crai
    touch ${sample}.g.vcf.gz
    touch ${sample}.g.vcf.gz.tbi
    touch ${sample}.duplicate.metrics.txt
    touch versions.yml
    """
}
