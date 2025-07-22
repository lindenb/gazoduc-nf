process RUN_ORAD {
label 'process_medium'
tag "${meta.id} ${ora_file.name}"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(oradir)
        tuple val(meta ),path(ora_file)
output:
        tuple val(meta),path("*q.gz"),emit:fastq
        path("versions.yml"),emit:versions
when:
        task.ext.when == null || task.ext.when
script:
        def sample = meta.id
"""
mkdir -p TMP

        ${oradir}/orad \
            --threads ${task.cpus} \\
            --path TMP \\
            --ora-reference "${oradir}/oradata" \\
            --name ${ora_file} \\
            -q

mv -v TMP/*.gz ./


    
cat <<-END_VERSIONS > versions.yml    
"${task.process}":
    orad: \$(${oradir}/orad --version | awk '(NR==1) {print \$3;}')
END_VERSIONS
"""

stub:
"""
touch ${meta.id}.R1.fq.gz
touch ${meta.id}.R2.fq.gz
touch versions.yml
"""
}
