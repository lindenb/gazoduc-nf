
process SMOOVE_ANNOTATE {
tag "${meta.id}"
label "process_single"
afterScript  "rm -rf TMP TMP2"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta4),path(gff3),path(gff3tbi)
	tuple val(meta ),path(vcf),path(vcfidx)
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def sample = task.ext.id?:meta.id
	def prefix = task.ext.prefix?:meta.id+".annot"
"""
	hostname 1>&2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP/TMP2
	export TMPDIR=\${PWD}/TMP/TMP2

    smoove annotate --gff ${gff3} "${vcf}" |\\
		bcftools sort \\
			--max-mem ${task.memory.giga}G  \\
			-T TMP/tmp \\
			-o "${prefix}.vcf.gz" \\
			-O z
	
	bcftools index --threads ${task.cpus} -f -t "${prefix}.vcf.gz"

cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}