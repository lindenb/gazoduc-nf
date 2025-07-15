
/**
 *
 * environ 20-30 minutes pour un WGS cram
 *
 */
process SMOOVE_CALL {
tag "${meta.id}"
label "process_single"
afterScript  "rm -rf TMP TMP2"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(exclude_bed)
	tuple val(meta ),path(bam),path(bai)
output:
	tuple path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
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
		--outdir TMP \\
		--exclude ${exclude_bed} \\
		--name ${sample} \\
		--fasta ${fasta} \\
		-p ${task.cpus} \\
		--genotype \\
		${bam}

    mv TMP/${sample}-smoove.genotyped.vcf.gz ./
    bcftools index -f -t ${sample}-smoove.genotyped.vcf.gz

    rm -rf TMP

cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}