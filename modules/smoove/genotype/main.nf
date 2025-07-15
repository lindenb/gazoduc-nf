process GENOTYPE_BAM {
tag "${row.sample}/${file(row.bam).name}/${file(merged).name}"
cache "lenient"
errorStrategy "retry"
maxRetries 5
memory "10g"
afterScript  "rm -rf TMP TMP2"
cpus 1 /* can only parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread. */
input:
	val(meta)
	val(genomeId)
	val(img)
	val(merged)
	val(row)
output:
	path("${row.sample}-smoove.regenotyped.vcf.gz"),emit:vcf
	path("${row.sample}-smoove.regenotyped.vcf.gz.tbi"),emit:tbi
	path("version.xml"),emit:version
script:
	def reference = params.genomes[genomeId].fasta
	def sample = row.sample
	def bam = row.bam
	def ref = file(reference)
	def vcf0 = file(merged)
	def xbam = file(bam)
"""
	hostname 1>&2
	#module load singularity/2.4.5
	module load bcftools/0.0.0
	mkdir -p TMP
	# smoove will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space
	mkdir -p TMP2
	export TMPDIR=\${PWD}/TMP2

	singularity exec\
		--home \${PWD} \
		--bind ${xbam.getParent()}:/bamdir \
		--bind ${vcf0.getParent()}:/mergeddir \
		--bind ${ref.getParent()}:/ref \
		--bind \${PWD}/TMP:/outdir \
		${img} \
		smoove genotype -x  \
			--vcf /mergeddir/${vcf0.name} \
			--outdir /outdir \
			--name ${sample} \
			--fasta /ref/${ref.name} \
			-p ${task.cpus} \
			/bamdir/${xbam.name}

	bcftools sort  --max-mem ${task.memory.giga}G  -T TMP -o ${sample}-smoove.regenotyped.vcf.gz -O z TMP/${sample}-smoove.genotyped.vcf.gz
	bcftools index -t -f ${sample}-smoove.regenotyped.vcf.gz


cat <<- EOF > versions.yml
"${task.process}":
    smoove: todo
EOF
"""
}