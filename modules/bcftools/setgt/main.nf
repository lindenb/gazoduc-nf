

process BCTOOLS_SETGT {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta),val(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
        path("versions.yml")
    script:
        def prefix = task.ext.prefix?:"${meta.id}.setgt"
        def args1 = task.ext.args1?:""
        if(args1.trim().isEmpty()) throw new IllegalArgumentException("do task.ext.args1 defined for ${task}")
    """
    mkdir -p TMP

    bcftools +setGT \\
        -O z \\
        -o "TMP/jeter.vcf.gz" \\
        "${vcf}" \\
        -- \\
        ${args1}

    bcftools index \\
        -f -t \\
        --threads "${task.cpus}" \\
        "TMP/jeter.vcf.gz"

    mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
    mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """

stub:
	def prefix = task.ext.prefix?:"${meta.id}.setgt"
        def args1 = task.ext.args1?:""
        if(args1.trim().isEmpty()) throw new IllegalArgumentException("do task.ext.args1 defined for ${task}")
"""
touch versions.yml ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi
"""
    }
