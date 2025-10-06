process BWA_INDEX {
tag "${fasta.name}"
label "process_medium"
conda "${moduleDir}/../../../conda/bwa.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(fasta)
output:
	tuple val(meta),path("${fasta.baseName}BWAIndex"),emit:bwa_index
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args = task.ext.args?:""
"""
    mkdir -p TMP
    bwa index ${args} -p TMP/${prefix} ${fasta}
    mv TMP "${fasta.baseName}BWAIndex"

cat << END_VERSIONS > versions.yml
${task.process}:
    bwa : todo
END_VERSIONS
"""

stub:
"""
mkdir -p TMP
touch TMP/${fasta.baseName}.amb
touch TMP/${fasta.baseName}.ann
touch TMP/${fasta.baseName}.bwt
touch TMP/${fasta.baseName}.pac
touch TMP/${fasta.baseName}.sa
mv TMP "${fasta.baseName}BWAIndex"
touch versions.yml
"""
}
