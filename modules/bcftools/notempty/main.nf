

process NOT_EMPTY_VCF {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    array 100
    when:
        task.ext.when == null || task.ext.when
    input:
        tuple val(meta),val(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*{.vcf.gz,.bcf}"),path("*{.tbi,.csi}"),optional:true,emit:vcf
        path("versions.yml")
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".notempty"
    """

    if test \$(bcftools index  -s ${vcf} | wc -l) -gt 0
    then 
        ln -s ${vcf} ${prefix}.${vcf.extension}
        ln -s ${vcfidx} ${prefix}.${vcfidx.extension}
    fi

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
END_VERSIONS
    """
    }