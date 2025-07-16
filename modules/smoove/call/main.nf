
/**
 *
 * environ 20-30 minutes pour un WGS cram
 *
 */
process SMOOVE_CALL {
tag "${meta.id}"
label "process_short"
afterScript  "rm -rf TMP TMP2"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta4),path(exclude_bed)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def sample = task.ext.id?:meta.id
"""
	hostname 1>&2
	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2


    smoove call \\
		--outdir TMP2 \\
		${exclude_bed?"--exclude \"${exclude_bed}\"":""} \\
		--name ${sample} \\
		--fasta ${fasta} \\
		-p ${task.cpus} \\
		--genotype \\
		${bam}


	find TMP2 1>&2

	bcftools index \\
		--threads ${task.cpus} \\
		-f -t \\
		TMP2/${sample}-smoove.genotyped.vcf.gz


    mv TMP2/${sample}-smoove.genotyped.vcf.gz ./
	mv TMP2/${sample}-smoove.genotyped.vcf.gz.tbi ./
    


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}