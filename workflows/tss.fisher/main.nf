workflow {
	TSS_FISHER(
		params.fasta,
		file(params.gtf),
		file(params.cases),
		file(params.controls)
		)
	}

workflow TSS_FISHER {
	take:
		fasta
		gtf
		cases
		controls
	main:

//				filter{T->T[0].matches("(chr)?[0-9]+")}.


		each_contig  = Channel.from(file(fasta+".fai")).
				splitCsv(sep:'\t',header:false).
				filter{T->T[0].matches("(chr)?[0-9]+")}.
				map{T->T[0]}
		each_contig.dump()
		
		sn_ch = SAMTOOLS_SN(Channel.of(
				[cases,"case"],
				[controls,"control"]
				))
		

		tss_bed_ch = TSS_CONTIG(gtf,each_contig)
		xtss_bed_ch = TSS_EXTEND_CONTIG(fasta, tss_bed_ch.output)


		ch1_ch = each_contig.combine(sn_ch.output.splitCsv(sep:'\t',header:false)).map{T->[
			contig: T[0],
			sample: T[1],
			bam : T[2],
			status: T[3]
			]}.combine(xtss_bed_ch.output.map{T->[
			contig: T[0],
			bed: T[1],
			extended: T[2]
			]}).
			filter{T->T[0].contig.equals(T[1].contig)}.
			map{T->T[0].plus(T[1])}.
			view()



		cov_ch = COVERAGE_BAM(fasta, ch1_ch)
		fisher_ch = FISHER(cov_ch.output.map{T->[T[0].contig,T[1]]}.groupTuple())
		MERGE(fisher_ch.output.collect())
	}

process SAMTOOLS_SN {
tag "${bams}"
label "process_low"
input:
	tuple path(bams),val(pheno)
output:
	path("sn.${pheno}.tsv"),emit:output
"""
module load samtools
cat "${bams}" | samtools samples | awk '{printf("%s\t%s\t${pheno}\\n",\$1,\$2);}' > sn.${pheno}.tsv
"""
}


process TSS_CONTIG {
tag "${contig}"
label "medium"
input:
	path(gtf)
	val(contig)
output:
	tuple val(contig),path("${contig}.tss.bed"),emit:output
script:
"""
${gtf.name.endsWith(".gz")?"gunzip -c":"cat"} '${gtf}' |\
	awk -F '\t' '(\$3=="gene" && \$1=="${contig}") {P=int(\$7=="+"?\$4:\$5);split(\$9,a,/["]/);printf("%s\t%d\t%d\t%s\t%d\\n",\$1,P-1,P,a[2],P-1);}' |\
	sort -T . -k1,1 -k2,2n > "${contig}.tss.bed"
"""
}

process TSS_EXTEND_CONTIG {
tag "${contig}"
label "medium"
input:
	val(fasta)
	tuple val(contig),path(bed)
output:
	tuple val(contig),path(bed),path("${contig}.tss.extend.bed"),emit:output
script:
	def binsize=500
"""
hostname 1>&2
module load bedtools

${bed.name.endsWith(".gz")?"gunzip -c":"cat"} '${bed}' |\
	bedtools slop -g '${fasta}.fai' -b {binsize} -i - |\
	sort -t '\t' -T . -k1,1 -k2,2n > ${contig}.tss.extend.bed
"""

}


process COVERAGE_BAM {
tag "${row.sample} ${row.contig}"
label "medium"
input:
	val(fasta)
	val(row)
output:
	tuple val(row),path("${row.sample}.${row.contig}.tsv"),emit:output
script:
"""
hostname 1>&2
module load samtools bedtools
set -o pipefail
mkdir -p TMP
set -x

samtools coverage --no-header --region "${row.contig}" "${row.bam}" |\\
	awk -F '\t' '{printf("%s\t0\t%d\t%s\\n",\$1,\$3,\$7);}' |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n > TMP/coverage.bed

head TMP/coverage.bed 1>&2


head '${row.extended}' 1>&2

samtools depth -r '${row.contig}' -a -b '${row.extended}'  --reference '${fasta}' '${row.bam}' |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	awk -F '\t' '{printf("%s\t%d\t%d\t%d\\n",\$1,int(\$2)-1,\$2,\$3);}' > TMP/depth.bed


head ${row.extended} TMP/depth.bed 1>&2


bedtools intersect -wa -wb -a '${row.extended}' -b TMP/depth.bed |\\
	sort -t '\t' -T TMP -k1,1 -k2,2n |\\
	bedtools intersect -wa -wb -a - -b TMP/coverage.bed |\\
		awk -F '\t' '{printf("%s\t%s\t%s\t%s\t%s\t%d\t%f\\n",\$1,\$2,\$3,\$4,\$5,int(\$7)-int(\$5),\$9 / \$13);}' > TMP/joined.bed


awk -F '\t' '(\$6==0) {SCORE=\$7;if(SCORE < 5) { S="NO";} else if(SCORE>10) {S="YES";} else {next;} printf("%s\t${row.status}\t%s\\n",\$4,S);}' TMP/joined.bed |\\
		sort -T TMP -t '\t' -k1,1 > "${row.sample}.${row.contig}.tsv"

"""
}


process FISHER {
tag "${contig} : ${L.size()}"
label "medium"
input:
	tuple val(contig), val(L)
output:
	path("output.${contig}.tsv"),emit:output
script:
"""
hostname
module load  R-project
mkdir -p TMP


cat ${L.join(" ")} > TMP/jeter.tsv

cat << '__EOF__' > TMP/jeter.R

data <- read.table("TMP/jeter.tsv", header = FALSE, sep = "\t", col.names = c("GENE", "STATUS", "CATEGORY"))
head(data)
contingency_table <- table(data\$GENE, data\$STATUS, data\$CATEGORY)
head(contingency_table)

result_df <- data.frame(
  GENE = rownames(contingency_table),
  CASE_YES = numeric(nrow(contingency_table)),
  CASE_NO = numeric(nrow(contingency_table)),
  CONTROL_YES = numeric(nrow(contingency_table)),
  CONTROL_NO = numeric(nrow(contingency_table)),
  fisher_test = numeric(nrow(contingency_table))
)

for (i in seq_along(result_df\$GENE)) {
  result_df[i, c("CASE_YES", "CASE_NO", "CONTROL_YES", "CONTROL_NO")] <- contingency_table[i, , ]
  
  fisher_result <- fisher.test(matrix(unlist(result_df[i, c("CASE_YES", "CASE_NO", "CONTROL_YES", "CONTROL_NO")]), nrow = 2, byrow = TRUE))
	
  result_df[i, "fisher_test"] <- fisher_result\$p.value

}

result_df <- result_df[result_df\$fisher_test <= 0.2, ]

result_df <- result_df[order(result_df\$fisher_test), ]

write.table( result_df , "output.${contig}.tsv", row.names=FALSE,quote=FALSE,sep='\t')

__EOF__

R --vanilla < TMP/jeter.R

"""
}

process MERGE {
executor "local"
input:
	val(L)
output:
	path("output.tsv"),emit:output
script:
"""
cat ${L.join(" ")} | LC_ALL=C sort -T . -t '\t' -k6,6g > output.tsv
"""
}
