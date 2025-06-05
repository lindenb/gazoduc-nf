/**
ORIGINAL snakemake workflow by Raphael Blanchet PhD.


*/

workflow {
	samplesheet_ch =Channel.fromPath(params.samplesheet).
		splitCsv(sep:'\t',header:true).
		map{[ [:], it.sample, it.R1, it.R2]}

	
	FASTQ_TO_BAM(
		file(params.fasta),
		file(params.fai),
		file(params.known_vcf),
		file(params.known_vcf_tbi),
		samplesheet_ch
		)
	}

process FASTQ_TO_BAM {
  tag "${sample} ${R1.name} ${R2.name}"
  // afterScript "rm -rf TMP"
  input:
	path(fasta)
	path(fai)
	path(known_vcf)
	path(known_vcf_tbi)
	tuple val(meta),val(sample),path(R1),path(R2)
  output:
	tuple val(meta),
		path("${sample}.cram"),
		path("${sample}.cram.crai"),
		path("${sample}.g.vcf.gz"),
		path("${sample}.g.vcf.gz.tbi"),
		emit:output
  script:
    def lib= meta.LIB?:sample
    def bwa="-M"
    def cpu_per_gpu = task.ext.cpu_per_gpu
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
                --gpusort \\
                --fix-mate  \\
                --knownSites "${known_vcf}" \\
                --out-recal-file "${sample}.bqsr.report.txt"  \\
		--gvcf \\
		--out-duplicate-metrics "${sample}.duplicate.metrics.txt" \\
                --tmp-dir TMP \\
                --logfile ${sample}.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(grep Parabr -m1 -A1 "${sample}.log" | grep -o 'Version [^ ]*' )
    END_VERSIONS
    """
}
