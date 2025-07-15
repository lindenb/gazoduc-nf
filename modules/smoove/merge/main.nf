
process SMOOVE_MERGE {
label "process_single"
tag "N=${meta.id}"
afterScript  "rm -rf TMP TMP2"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path("VCFS/*")
output:
	tuple val(meta),path("*.sites.vcf.gz"),path("*.sites.vcf.gz.tbi"),emit:tbi
	path("versions.yml"),emit:versions

script:
	def prefix = task.ext.prefix?:meta.id+".sites"
"""
	hostname 1>&2

	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space


	smoove merge \
		--name merged \
		--outdir TMP \
		--name merged \
		--fasta ${fasta} \
		`find VCFS -name "*.vcf.gz" |sort -T TMP`


	bcftools sort  \\
		--max-mem ${task.memory.giga}G  \\
		-T TMP/sort \\
		-o ${prefix}.vcf.gz \\
		-O z \\
		TMP/merged.sites.vcf.gz
	
	bcftools index -t -f ${prefix}.vcf.gz


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}
