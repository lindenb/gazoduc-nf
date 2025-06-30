process GATK_POSSIBLE_DENOVO {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
        tuple val(meta4),path(pedigree)
       tuple val(meta),path(vcf),path(vcfidx)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def prefix = task.ext.prefix?:vcf.baseName+".denovo"
        def input_is_bcf = vcf.name.endsWith(".bcf")
    """
    mkdir -p TMP

    if ${input_is_bcf}
    then
        bcftools view -O z --threads ${task.cpus} -o TMP/jeter1.vcf.gz '${vcf}'
        bcftools index -f -t --threads ${task.cpus} TMP/jeter1.vcf.gz
    fi


    awk '{SEX=\$5;if(SEX=="male") SEX="1"; if(SEX=="female") SEX="2"; PHENO="0";
         if(\$6=="affected" || \$6=="case") PHENO="2";
         if(\$6=="unaffected" || \$6=="control") PHENO="1";
         printf("%s\t%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$4,SEX,PHENO)}' "${pedigree}" > TMP/jeter.ped

    gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantAnnotator \\
        -R "${fasta}" \\
        --annotation PossibleDeNovo \\
        --pedigree  TMP/jeter.ped \\
        -V ${input_is_bcf? "TMP/jeter1.vcf.gz" : "\"${vcf}\""} \\
        -O "TMP/${prefix}.vcf.gz"

    rm -f TMP/jeter1.vcf.gz TMP/jeter1.vcf.gz.tbi

    bcftools index -f -t \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}.vcf.gz"
    
    mv "TMP/${prefix}.vcf.gz" ./
    mv "TMP/${prefix}.vcf.gz.tbi" ./
    """
    }
