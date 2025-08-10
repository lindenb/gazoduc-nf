process DUPHOLD {
tag "${meta.id?:""}"
label 'process_short'
conda "${moduleDir}/../../conda/duphold.yml"
afterScript "rm -rf TMP"

input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta), path(bam), path(bai), path(svvcf), path(optional_snp), path(optional_snp_idx)
output:
    tuple val(meta), path("*.vcf.gz") ,path("*.tbi"),emit: vcf
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
"""
hostname 1>&2
mkdir -p TMP

duphold \\
    ${args} \\
    ${optional_snp ? "--snp \"${optional_snp}\"" : ""} \\
    --threads ${task.cpus} \\
    --output TMP/jeter.vcf.gz \\
    --vcf "${svvcf}" \\
    --bam "${alignment_file}" \\
    --fasta "${fasta}"

bcftools index -f -t --threads "${task.cpus}" TMP/jeter.vcf.gz

mv  TMP/jeter.vcf.gz ${prefix}.vcf.gz
mv  TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    duphold: \$(duphold -h | head -n 1 | sed -e "s/^version: //")
END_VERSIONS
"""
}
