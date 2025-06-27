include {k1_signature} from '../../../modules/utils/k1.nf'

workflow TRIOS {
    take:
        meta
        fasta
        fai
        dict
        vcf
        pedigree
    main:
      BCTOOLS_MENDELIAN(fai, vcf, pedigree)
      GATK_POSSIBLE_DENOVO(fasta,fai,dict,BCTOOLS_MENDELIAN.out.vcf,pedigree)
}


process BCTOOLS_MENDELIAN {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"

    input:
        tuple val(meta1),path(fai)
        tuple val(meta),val(contig),path(vcfidx)
        tuple val(meta2),path(pedigree)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def prefix = meta.id+".md."+contig
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

    bcftools +mendelian2 \\
        --threads  "${task.cpus}" \\
        -O z \\
        -o "TMP/${prefix}.vcf.gz" \\
        `cat TMP/jeter.rule` \\
        -m a -P '${pedigree}' "${vcf}"

    bcftools index \\
        -ft \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}.vcf.gz"

    mv "TMP/${prefix}.vcf.gz" ./
    mv "TMP/${prefix}.vcf.gz.tbi" ./
    """
    }


process GATK_POSSIBLE_DENOVO {
    label "process_single"
    tag "${meta.id}"
    conda "${moduleDir}/../../../conda/bioinfo.01.yml"

    input:
        tuple val(meta1),path(fasta)
        tuple val(meta2),path(fai)
        tuple val(meta3),path(dict)
        tuple val(meta),path(vcf),path(vcfidx)
        tuple val(meta4),path(pedigree)
    output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    script:
        def prefix = task.ext.prefix?:vcf.simpleName+".denovo"
    """
    mkdir -p TMP

    gatk --java-options "-XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantAnnotator \\
        -R "${fasta}" \\
        --annotation PossibleDeNovo \\
        --pedigree  '${pedigree}' \\
        -V "${vcf}" \\
        -O "TMP/${prefix}.vcf.gz"

    bcftools index -f -t \\
        --threads "${task.cpus}" \\
        "TMP/${prefix}.vcf.gz"
    
    mv "TMP/${prefix}.vcf.gz" ./
    mv "TMP/${prefix}.vcf.gz.tbi" ./
    """
    }

