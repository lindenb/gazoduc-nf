workflow VEP {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        if(isGRCH38(fai)) {
            VEP_GRCH38(meta,fasta,fai,dict,vcf)
            vcf = VEP_GRCH38.out.vcf
        } else {
            throw new IllegalArgumentException("VEP: unknown build");
        }
}

workflow VEP_GRCH38 {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        
    emit:
        vcf
}


process ANNOT_VEP_GRCH38 {
tag "${sample} ${vcf.baseName}"
label "process_quick"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
        path(fasta)
        path(fai)
        path(fix_vep)
        tuple val(sample),path(vcf),path(tbi)
output:
        tuple val(sample),path("${vcf.baseName}.vep.vcf.gz"),path("${vcf.baseName}.vep.vcf.gz.tbi"),emit:output
script:
"""
mkdir -p TMP

bcftools view "${vcf}" |\\
vep \\
    --cache --dir "${params.vep}" \\
    --species homo_sapiens \\
    --merged \\
    --assembly GRCh38 \\
    --cache_version "${params.vep_cache_version}" \\
    --fasta "${params.fasta}" \\
    --offline \\
    --symbol \\
    --format vcf \\
    --force_overwrite \\
    --sift=b \\
    --polyphen=b \\
    --custom "${params.gnomad},gnomADg,vcf,exact,0,AC,AC_afr,AC_afr_female,AC_afr_male,AC_ami,AC_ami_female,AC_ami_male,AC_amr,AC_amr_female,AC_amr_male,AC_asj,AC_asj_female,AC_asj_male,AC_eas,AC_eas_female,AC_eas_male,AC_female,AC_fin,AC_fin_female,AC_fin_male,AC_male,AC_nfe,AC_nfe_female,AC_nfe_male,AC_oth,AC_oth_female,AC_oth_male,AC_raw,AC_sas,AC_sas_female,AC_sas_male,AF,AF_afr,AF_afr_female,AF_afr_male,AF_ami,AF_ami_female,AF_ami_male,AF_amr,AF_amr_female,AF_amr_male,AF_asj,AF_asj_female,AF_asj_male,AF_eas,AF_eas_female,AF_eas_male,AF_female,AF_fin,AF_fin_female,AF_fin_male,AF_male,AF_nfe,AF_nfe_female,AF_nfe_male,AF_oth,AF_oth_female,AF_oth_male,AF_raw,AF_sas,AF_sas_female,AF_sas_male,nhomalt,nhomalt_afr,nhomalt_afr_female,nhomalt_afr_male,nhomalt_ami,nhomalt_ami_female,nhomalt_ami_male,nhomalt_amr,nhomalt_amr_female,nhomalt_amr_male,nhomalt_asj,nhomalt_asj_female,nhomalt_asj_male,nhomalt_eas,nhomalt_eas_female,nhomalt_eas_male,nhomalt_female,nhomalt_fin,nhomalt_fin_female,nhomalt_fin_male,nhomalt_male,nhomalt_nfe,nhomalt_nfe_female,nhomalt_nfe_male,nhomalt_oth,nhomalt_oth_female,nhomalt_oth_male,nhomalt_raw,nhomalt_sas,nhomalt_sas_female,nhomalt_sas_male" \\
    -o TMP/jeter.vcf \\
    --quiet \\
    --vcf \\
    --no_stats \\
    --plugin SpliceAI,snv=${params.vep_spliceai_snv},indel=${params.vep_spliceai_indel} \\
    --plugin LOEUF,file=${params.vep_loeuf},match_by=gene 

bcftools view -O z -o "${vcf.baseName}.vep.vcf.gz" TMP/jeter.vcf
bcftools index -t -f "${vcf.baseName}.vep.vcf.gz"
"""
}
