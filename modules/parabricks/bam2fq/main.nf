process PB_BAM2FQ {
  tag "${meta.id}"
  label 'process_short'

  afterScript "rm -rf TMP"
  input:
    tuple val(meta),path(bam),path(fasta),path(fai)
  output:
      tuple val(meta), path("*.gz")  , emit: fastqs
      path("${meta.id}.bam2fq.log")
      path("versions.yml"),emit:versions
  when:
        task.ext.when == null || task.ext.when
  script:
    def sample = meta.id?bam.simpleName
    def with_orphan = task.ext.with_orphan?:false
    def with_singleton = task.ext.with_singleton?:false
    def cpu_per_gpu = task.ext.cpu_per_gpu
    def bqsr_args = with_bqsr && known_vcf? "--knownSites \"${known_vcf.name}\"   --out-recal-file \"${sample}.bqsr.report.txt\"   " : ""
    def low_memory = task.ext.low_memory==true || task.attempt>1 ? "--low-memory" : ""
    def fixmate_args = (task.ext.with_fixmate==true?"--fix-mate":"")
 """
	mkdir -p TMP/OUT

  pbrun bam2fq \\
      --num-gpus ${task.ext.gpus} \\
      --num-threads ${task.cpus} \\
      --num-cpu-threads-per-stage ${cpu_per_gpu} \\
      --memory-limit ${task.memory.giga} \\
      --ref ${fasta} \\
      --out-prefix "TMP/OUT/${sample}." \\
      --out-suffixF  .R1.fastq.gz \\
      --out-suffixF2 .R2.fastq.gz \\
      ${with_orphan   ?"--out-suffixO  .orphan1.fastq.gz":""} \\
      ${with_orphan   ?"--out-suffixO2 .orphan2.fastq.gz":""} \\
      ${with_singleton?"--out-suffixS  .singletons.fastq.gz":""} \\
      --in-bam "${bam}" \\
      --tmp-dir TMP \\
      --logfile ${sample}.bam2fq.log \\
      ${low_memory}

find .  1>&2

mv "TMP/OUT/*.gz ./


cat << END_VERSIONS > versions.yml
"${task.process}":
    pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
END_VERSIONS
"""


stub:
    def sample = meta.id
    """
    touch ${sample}.R1.fastq.gz
    touch ${sample}.R2.fastq.gz
    touch ${sample}.bam2fq.txt
    touch versions.yml
    """
}
