process SCATTER_INTERVALS_BY_NS {
tag "${meta1.id?:fasta.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
    tuple val(meta1),path("*.interval_list"),emit:interval_list
    path("versions.yml"),emit:versions
script:
    def max_to_merge = task.ext.max_to_merge?:""
    def output_type = task.ext.output_type?:""
    if(max_to_merge.toString().trim().isEmpty()) throw new IllegalArgumentException("undefined max_to_merge for ${task.process}");
    if(output_type.toString().trim().isEmpty()) throw new IllegalArgumentException("undefined output_type for ${task.process}");
    def prefix = task.ext.prefix?:fasta.baseName
"""
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ScatterIntervalsByNs \\
    --REFERENCE "${fasta}" \\
    --MAX_TO_MERGE ${max_to_merge} \\
    --OUTPUT "TMP/jeter.interval_list" \\
    --OUTPUT_TYPE ${output_type}

mv "TMP/jeter.interval_list" ./${prefix}.interval_list

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
"""
}
