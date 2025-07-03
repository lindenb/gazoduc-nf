
process SEQTK_SPLIT {
tag "${meta.id?:R1.name}"
afterScript "rm -rf TMP OUT/fastq.tsv"
// cpus 4 set in config file
input:
	val(meta),path(R1),path(R2)
	val(row)
output:
	tuple val(meta),path("fastq.tsv"),emit:output
	path("versions.yml"),emit:versions
script:
	def COL = row.grep({!(it.key.equals("R1") || it.key.equals("R2"))}).collect{it.key}.join("\t")
	def ROW = row.grep({!(it.key.equals("R1") || it.key.equals("R2"))}).collect{it.value}.join("\t")

	seqtk = params.seqtk.split2.args

if(row.containsKey("R2") && !row.R2.isEmpty() && !row.R2.equals("."))
"""
hostname 1>&2
${moduleLoad("seqkit")}
set -o pipefail

mkdir -p TMP

seqkit split2 -O TMP --force \\
    -j ${task.cpus} ${seqtk} \\
    -1 '${row.R1}' \\
    -2 '${row.R2}'
mv -v TMP OUT


find \${PWD}/OUT/ -type f -name "*q.gz" |\
	awk -F '.' '{printf("%s\t%s\\n",\$(NF-2),\$0);}' |\
	sort -T . -t '\t' -k1,1V |\
	paste - - |\
	cut -f 2,4 |\
	awk -F '\t' '{printf("${ROW}\t%s\\n",\$0);}' > OUT/fastq.tsv

test -s OUT/fastq.tsv

echo "${COL}\tR1\tR2" > fastq.tsv
cat OUT/fastq.tsv >> fastq.tsv




cat << EOF > versions.yml
${task.process}:
    seqkit: "\$(seqkit  version)"
EOF
"""
}