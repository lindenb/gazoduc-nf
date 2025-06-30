include {k1_signature} from '../../../modules/utils/k1.nf'


process BCTOOLS_MENDELIAN2 {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"
    afterScript "rm -rf TMP"
    input:
        tuple val(meta1),path(fai)
        tuple val(meta2),path(pedigree)
        tuple val(meta),val(contig),path(vcfidx)
    output:
        tuple val(meta),path("*{.vcf.gz,.bcf}"),path("*{.tbi,.csi}"),emit:vcf
    script:
        def prefix = meta.id+".mendelian."+contig
        def suffix = ".bcf"
    """
    mkdir -p TMP


    cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
    1:${k1.hg38}\t--rules GRCh38
    1:${k1.hg19}\t--rules GRCh37
    EOF

    awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
        sed 's/^chr//' |\\
        sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

    join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
        sort |\\
        uniq > TMP/jeter.rule


    awk '{SEX=1;if(\$5=="female")) SEX=2;printf("%s\t%s\t0\t0\t%s\\n",\$1,\$1,SEX)}' ${pedigree} >> TMP/jeter.ped


    bcftools +mendelian2 \\
        --threads  "${task.cpus}" \\
        -O ${suffix.contains("bcf")?"b":"z"} \\
        -o "TMP/${prefix}${suffix}" \\
        `cat TMP/jeter.rule` \\
        -m a -P TMP/jeter.ped "${vcf}"

    bcftools index \\
        -f \\
        ${suffix.contains("bcf")?"":"-t"} \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}}${suffix}"

    mv "TMP/${prefix}.*" ./
    """
    }