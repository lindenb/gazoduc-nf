
process SMOOVE_PASTE {
tag "${meta.id}"
label "process_single"
afterScript  "rm -rf TMP TMP2"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path("VCFS/*")
output:
	tuple val(meta ),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:meta.id+".paste"
"""
	hostname 1>&2
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2

	smoove paste --name "${prefix}" VCFS/*.vcf.gz

	bcftools index -t -f "${prefix}.smoove.square.vcf.gz"

cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}
