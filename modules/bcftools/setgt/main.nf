

process BCTOOLS_SETGT {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    when:
        task.ext.when == null || task.ext.when
    input:
        tuple val(meta),val(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*{.vcf.gz,.bcf}"),path("*{.tbi,.csi}"),emit:vcf
        path("versions.yml")
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".setgt"
        def suffix = task.ext.suffix?:".bcf"
        def args1 = task.ext.args1?:""
        if(args1.trim().isEmpty()) throw new IllegalArgumentException("do task.ext.args1 defined for ${task}")
    """
    mkdir -p TMP

    bcftools +setGT \\
        -O ${suffix.contains("bcf")?"b":"z"} \\
        -o "TMP/${prefix}${suffix}" \\
        "${vcf}" \\
        -- \\
        ${args1}

    bcftools index \\
        -f \\
        ${suffix.contains("bcf")?"":"-t"} \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}${suffix}"

    mv TMP/${prefix}* ./


cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """
    }