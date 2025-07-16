process SMOOVE_GENOTYPE {
label "process_short"
tag "${meta.id}"
afterScript  "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(sitesvcf),path(sitesvcfidx)
	tuple val(meta ),path(bam),path(bai)

output:
	tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def sample=meta.id
	def prefix = task.ext.prefix?:meta.id+".genotype"
"""
	hostname 1>&2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP/TMP2
	export TMPDIR=\${PWD}/TMP/TMP2

	smoove genotype -d -x \\
		--vcf "${sitesvcf}" \\
		--outdir TMP \\
		--name ${sample} \\
		--fasta ${fasta} \\
		-p ${task.cpus} \\
		${bam}

	bcftools sort  \\
		--max-mem ${task.memory.giga}G \\
		-T TMP/sort \\
		-o ${prefix}.vcf.gz \\
		-O z \\
		TMP/${sample}-smoove.genotyped.vcf.gz
	
	bcftools index --threads ${task.cpus} -t -f ${prefix}.vcf.gz


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}