
process SMOOVE_MERGE {
label "process_short"
tag "${meta.id}"
afterScript  "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta ),path("VCFS/*")
output:
	tuple val(meta),path("*.sites.vcf.gz"),path("*.sites.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:meta.id+".sites"
"""
	hostname 1>&2

	mkdir -p TMP/TMP2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	export TMPDIR=\${PWD}/TMP/TMP2


	smoove merge \\
		--outdir TMP \\
		--name merged \\
		--fasta ${fasta} \\
		`find VCFS -name "*.vcf.gz" |sort -T TMP`


	bcftools sort  \\
		--max-mem ${task.memory.giga}G  \\
		-T TMP/sort \\
		-o ${prefix}.vcf.gz \\
		-O z \\
		TMP/merged.sites.vcf.gz
	
	bcftools index \\
		--threads ${task.cpus} \\
		-t -f ${prefix}.vcf.gz


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}
