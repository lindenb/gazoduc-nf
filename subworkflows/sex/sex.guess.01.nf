include {getKeyValue } from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES01 } from '../../modules/samtools/samtools.samples.01.nf'

workflow SEX_GUESS_01 {
	take:
		meta
		reference
		bams
	main:
		xy_ch = sexualContigs(meta,reference).out.splitCsv(header:false,sep:'\t')
		sn_bam_ch = SAMTOOLS_SAMPLES01(meta, reference, bams).splitCsv(header:false,sep:'\t')
		sex_count_ch = sexualChromosomeCount(meta, reference, sn_bam_ch.combine(xy_ch))
		sn_sex_ch  = sampleSex(meta, sex_count_ch.collect())
	emit:
		out = sn_sex_ch.out
		version = sn_sex_ch.version
	}

process sexualContigs {
executor "local"
input:
	val(meta)
	val(reference)
output:
	path("XY.bed"),emit:out
script:
"""
test -s "${reference}.fai"
awk -F '\t' '(\$1 ~ /^(chr)?X\$/) {printf("%s\t%s\tX\\n",\$1,\$2);}' "${reference}.fai" > X.bed
test -s X.tsv
awk -F '\t' '(\$1 ~ /^(chr)?Y\$/) {printf("%s\t%s\tY\\n",\$1,\$2);}' "${reference}.fai" > Y.bed
test -s Y.tsv

cat X.tsv Y.tsv > XY.tsv
rm X.tsv Y.tsv

test `wc -l < XY.tsv` -eq 2
"""
}

process sexualChromosomeCount {
tag "${sample} ${contig} ${file(bam).name}"
label "process_low"
input:
	val(meta)
	val(reference)
	tuple val(sample),val(bam),val(contig),val(contigLen),val(contigType)
output:
	path("count.txt"),emit:count
script:
	def mapq = getKeyValue(meta,"mapq","30")
"""
hostname 2>&1
module load ${getModules("samtools")}

samtools view -@ ${task.cpus} -c  -F 3844 --min-MQ ${mapq} \
	--reference "${reference}" \
	"${bam}" "${contig}" |\
	awk '{printf("${sample}\t${bam}\t${reference}\t${contigType}\t${contigLen}\t%s\\n",\$1);}' > count.txt
"""
}

process  sampleSex {
executor "local"
tag "${KEY[0]} ${file(KEY[1]).name}"
input:
	val(meta)
	tuple val(KEY),val(L)
output:
	path("samples.sex.tsv"),emit:out
	path("version.xml"),emit:version
script:
	def treshold = getKeyValue(meta,"treshold","10")
"""
cat << EOF > jeter.txt
${L.join("\n")}
EOF

xargs -a jeter.txt -L 1 cat | sort -T . -t '\t' -k1,1 -k4,4 |\
	awk -F '\t' '(\$4=="X") {FX=int(\$6)/(int(\$5)*1.0);next;} (\$4=="Y") {S="male";FY=int(\$6)/(int(\$5)*1.0);if(FX > FY * (${treshold}*1.0)) {S="female"};printf("%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,S);next;}' > samples.sex.tsv

rm jeter.txt

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">extract count from sexual chromosome, guess the sex</entry>
	<entry key="treshold">${treshold}</entry>
</properties>
EOF
"""
}
