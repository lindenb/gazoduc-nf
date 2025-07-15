
/**
 *
 * environ 20-30 minutes pour un WGS cram
 *
 */
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


    smoove annotate \\
		--outdir TMP \\
		--exclude ${exclude_bed} \\
		--name ${sample} \\
		--fasta ${fasta} \\
		-p ${task.cpus} \\
		--genotype \\
		${bam}




	# update dict
	java -jar -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP \${PICARD_JAR} UpdateVcfSequenceDictionary \
		I=TMP/jeter1.vcf.gz O=TMP/jeter2.vcf.gz SD=${reference}
	mv TMP/jeter2.vcf.gz TMP/jeter1.vcf.gz

	# sort vcf
	bcftools sort  --max-mem ${task.memory.giga}G  -T TMP/tmp -o "${prefix}smoove.bcf" -O b  TMP/jeter1.vcf.gz
	bcftools index -f "${prefix}smoove.bcf"


    rm -rf TMP

cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}