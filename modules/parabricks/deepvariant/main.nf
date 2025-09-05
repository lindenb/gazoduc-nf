
process PB_DEEPVARIANT {
  tag "${meta.id}"
  label 'process_short'

  afterScript "rm -rf TMP"
  input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta),path(bam),path(bai)
  output:
	tuple val(meta), path("${meta.id}.g.vcf.gz"),path("${meta.id}.g.vcf.gz.tbi"),emit:gvcf, optional:true
	path("${meta.id}.deepvariant.log")
	path("versions.yml"),emit:versions
  script:
    def sample = meta.id
    def args1 = task.ext.args1?:""
 """
	mkdir -p TMP

  pbrun deepvariant \\
      --num-gpus ${task.ext.gpus} \\
      --ref ${fasta} \\
      --in-bam "${bam}" \\
      --gvcf \\
      --out-variants "${sample}.g.vcf.gz" \\
      --tmp-dir TMP \\
      --logfile ${sample}.deepvariant.log \\
      ${args1}


cat << END_VERSIONS > versions.yml
"${task.process}":
    pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
END_VERSIONS
"""


stub:
def sample = meta.id
"""
  touch ${sample}.g.vcf.gz
  touch ${sample}.g.vcf.gz.tbi
  touch ${sample}.log
  touch ${sample}.deepvariant.log
  touch versions.yml
"""
}
