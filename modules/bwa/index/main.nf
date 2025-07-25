process BWA_INDEX {
tag "${fasta.name}"
label "process_medium"
conda "${moduleDir}/../../../conda/bwa.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(fasta)
output:
	tuple val(meta),path("${fasta.baseName}BWAIndex"),emit:bwa_index
    path("versions.yml")
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
}
