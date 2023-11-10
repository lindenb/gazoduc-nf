include {moduleLoad} from '../../modules/utils/functions.nf'
include {GET_LIFTOVER_CHAINS_02} from '../../modules/ucsc/liftover.chains.02.nf'

workflow LIFTOVER_BAM {
	take:
		meta
		genomeIdDest
		row // contains sample,bam
	main:
		version_ch = Channel.empty()


		chain_ch = GET_LIFTOVER_CHAINS_02(
			[:],
			row.map{T->params.genomes[T.genomeId].ucsc_name}.unique(),
			params.genomes[genomeIdDest].ucsc_name
			)

		ch1 = row.combine(chain_ch.output).
			filter{T->params.genomes[T[0].genomeId].ucsc_name.equals(T[1])}.
			map{T->T[0].plus("chain":T[3])}



		contigs_ch = BAM2CONTIGS([:],ch1)

		ch2 = contigs_ch.output.splitCsv(header:true,sep:',').
			map{T->T[0].plus("contig":T[1].contig)}

		lift_ch = LIFT_BAM([:], genomeIdDest, ch2)


		merged_ch = MERGE_LIFTED([:],lift_ch.groupTuple())

		rows = merged_ch.output.map{T->[
			"sample":T[0],
			"genomeId":genomeIdDest,
			"bam":T[1]
			]}
	emit:
		rows = rows

}


	



process BAM2CONTIGS {
tag "${row.bam}"
afterScript "rm -rf TMP"
cpus 1
input:
      	val(meta)
        val(row)
output:
       	tuple val(row),path("contigs.csv"),emit:output
script:
       	def prefix = row.prefix?:""
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

samtools idxstats --threads ${task.cpus} "${row.bam}" | awk -F '\t' 'BEGIN {print "contig" } (\$1!="*") {printf("%s\\n",\$1);}' > contigs.csv
"""
}

process LIFT_BAM {
tag "sample:${row.sample} chrom:${row.contig}"
afterScript "rm -rf TMP"
memory "5g"
cpus 3
input:
	val(meta)
	val(genomeIdDest)
	val(row)
output:
	tuple val("${row.sample}"),path("contig.bam"),emit:output
when:
	!row.contig.equals("*")
script:
	def fasta1 = params.genomes[row.genomeId].fasta
	def fasta2 = params.genomes[genomeIdDest].fasta
	def smallerBam = true
"""
hostname 1>&2
mkdir -p TMP
${moduleLoad("samtools jvarkit")}
set -o pipefail
set -x

java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar convertliftoverchain \
	-R1 '${fasta1}' \
	-R2 '${fasta2}' \
	'${row.chain}' > jeter.chain

samtools view  -O BAM --exclude-flags 4 --uncompressed -T '${fasta1}' "${row.bam}" "${row.contig}" |\
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bamliftover \
		--chain jeter.chain \
		-R "${fasta1}" \
		-R2 "${fasta2}" \
		--drop-seq --unmapped \
		--samoutputformat BAM -o TMP/jeter.bam


samtools sort  --threads ${task.cpus} -T TMP/tmp -m '${task.memory.giga}G' -o TMP/jeter2.bam TMP/jeter.bam 

mv TMP/jeter2.bam contig.bam
rm jeter.chain
"""
}

process MERGE_LIFTED {
tag "${sample} N=${L.size()}"
input:
	val(meta)
	tuple val(sample),val(L)
output:
	tuple val(sample),path("${sample}.lift.bam"),path("${sample}.lift.bam.bai"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools")}

cat << EOF > jeter.list
${L.join("\n")}
EOF

samtools merge -o "${sample}.lift.bam" --write-index --threads ${task.cpus} -b jeter.list
rm jeter.list
"""
}
