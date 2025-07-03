process GATK_CALCULATE_GENOTYPE_POSTERIORS {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
        tuple val(meta4),path(pedigree)
        tuple val(meta5),path(optional_supporting_vcf),path(optional_supporting_vcf_tbi)
       tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
        path("versions.yml"),emit:versions
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".denovo"
        def input_is_bcf = vcf.name.endsWith(".bcf")
        def args1 = task.ext.args1?:""
        def supporting = optional_supporting_vcf?" -supporting ${optional_supporting_vcf}":""
    """
    mkdir -p TMP

    if ${input_is_bcf}
    then
        bcftools view -O z --threads ${task.cpus} -o TMP/jeter1.vcf.gz '${vcf}'
        bcftools index -f -t --threads ${task.cpus} TMP/jeter1.vcf.gz
    fi


    awk -f "${moduleDir}/../possibledenovo/pedigree4gatk.awk' "${pedigree}" > TMP/jeter.ped

    gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CalculateGenotypePosteriors \\
        -R "${fasta}" \\
        ${args1} \\
        ${supporting} \\
        --pedigree  TMP/jeter.ped \\
        -V ${input_is_bcf? "TMP/jeter1.vcf.gz" : "\"${vcf}\""} \\
        -O "TMP/${prefix}.vcf.gz"

    rm -f TMP/jeter1.vcf.gz TMP/jeter1.vcf.gz.tbi

    bcftools index -f -t \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}.vcf.gz"
    
    mv "TMP/${prefix}.vcf.gz" ./
    mv "TMP/${prefix}.vcf.gz.tbi" ./

cat << END_VERSIONS > versions.yml
"${task.process}":
	gatk: "\$(gatk --version 2>&1  | paste -s -d ' ' | tr -c -d 'A-Za-z0-9._-' )"
END_VERSIONS
    """
    }
