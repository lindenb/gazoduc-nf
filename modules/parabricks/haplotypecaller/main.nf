
process PB_HAPLOTYPECALLER {
  tag "${meta.id}"
  label 'process_short'

  afterScript "rm -rf TMP"
  input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path(bam),path(bai),path(recal_file)
  output:
	tuple val(meta),
		path("${meta.id}.g.vcf.gz"),
		path("${meta.id}.g.vcf.gz.tbi"),emit:bams
	path("${meta.id}.hc.log"),emit:log
	path("versions.yml"),emit:versions
  script:
    def sample = meta.id
    def recal_args = recal_file?"--in-recal-file \"${recal_file.name}\"":""
 """
	mkdir -p TMP

           pbrun haplotypecaller \\
                --num-gpus ${task.ext.gpus} \\
                --ref ${fasta} \\
                --in-bam "${bam}" \\
                --gvcf \\
                --out-variants "${sample}.g.vcf.gz" \\
                --tmp-dir TMP \\
                --logfile ${sample}.hc.log \\
		${recal_args}


cat <<-END_VERSIONS > versions.yml
"${task.process}":
        pbrun: \$(grep Parabr -m1 -A1 "${sample}.log" | grep -o 'Version [^ ]*' )
END_VERSIONS
"""

 stub:
    def sample = meta.id
    """
    touch ${sample}.g.vcf.gz
    touch ${sample}.g.vcf.gz.tbi
    touch ${sample}.hc.log
    touch versions.yml
    """
}
