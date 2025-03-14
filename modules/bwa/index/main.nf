process BWA_INDEX {
tag "${fasta.name}"
label "process_medium"
conda "${moduleDir}/../../../conda/bwa.yml"
afterScript "rm -rf TMP"
input:
	path(fasta)
output:
	path("BWA"),emit:output
script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
"""
    mkdir -p TMP
    bwa index -p TMP/${prefix} ${fasta}
    mv TMP BWA
"""
}
