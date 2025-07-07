include {VEP_INSTALL_PLUGINS} from '../../../modules/vep/install.plugins/main.nf'

workflow VEP {
    take:
        meta
        fasta
        fai
        dict
        vcf
    main:
        versions = Channel.empty()
        VEP_INSTALL_PLUGINS([[:],"SpliceAI,LOEUF,UTRAnnotator"])
        versions = versions.mix(VEP_INSTALL_PLUGINS.out.versions)

        DOWNLOAD_UTR_ANNOTATOR(fasta)
        versions = versions.mix(DOWNLOAD_UTR_ANNOTATOR.out.versions)
        
        APPLY_VEP(
            fasta,
            fai,
            VEP_INSTALL_PLUGINS.out.directory,
            DOWNLOAD_UTR_ANNOTATOR.out.output,
            vcf
            )
        versions = versions.mix(APPLY_VEP.out.versions)
    emit:
        vcf=APPLY_VEP.out.vcf
        versions
}

process DOWNLOAD_UTR_ANNOTATOR {
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
output:
    tuple val(meta1),path("*.txt"),emit:output
    path("versions.yml"),emit:versions
script:
    def f1="uORF_5UTR_GRCh38_PUBLIC.txt"
"""
wget -O "${f1}" "https://github.com/ImperialCardioGenetics/UTRannotator/raw/refs/heads/master/${f1}"


cat << END_VERSIONS > versions.yml
"${task.process}":
    db: todo
END_VERSIONS
"""
}

process APPLY_VEP {
tag "${meta.id?:vcf.baseName}"
label "process_single"// single is killed
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(plugindir)
    tuple val(meta4),path(utrfile)
    tuple val(meta),path(vcf),path(tbi)
output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
        path("versions.yml"),emit:versions
script:
    def assembly = task.ext.assembly?:"GRCh38"
    def species = task.ext.species?:"homo_sapiens"
    def merged =  task.ext.merged?:true
    def vep_directory = task.ext.vep_directory?:"/LAB-DATA/GLiCID/projects/BiRD_resources/apps/vep"
    def prefix = task.ext.prefix?:vcf.baseName+".vep"
    def cache_version = task.ext.cache_version?:"113"
    def spliceai_snv= task.ext.spliceai_snv?:"/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz"
    def spliceai_indel = task.ext.spliceai_index?:"/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz"
    def loeuf = task.ext.loeuf?:"/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/LOEUF/supplement/loeuf_dataset_grch38.tsv.gz"
    def gnomad = task.ext.gnomad?:"/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/broadinstitute.org/gnomad/3.0/gnomad.genomes.r3.0.sites.vcf.gz"
"""
mkdir -p TMP

bcftools query -l '${vcf}' |\\
    awk 'END{N=NR;if(N==0) N=1;B=5000.0/N;if(B<1000.0) B=1000;printf("%d\\n",int(B));}' > TMP/jeter.buffer_size


bcftools view "${vcf}" |\\
vep \\
    ${task.cpus>1?"--fork ${task.cpus}":""} \\
    --buffer_size `cat TMP/jeter.buffer_size` \\
    --cache \\
    --dir "${vep_directory}" \\
    --species ${species} \\
    ${merged?"--merged":""} \\
    --assembly ${assembly} \\
    --cache_version "${cache_version}" \\
    --fasta "${params.fasta}" \\
    --dir_plugins "${plugindir}/" \\
    --offline \\
    --symbol \\
    --format vcf \\
    --force_overwrite \\
    --sift=b \\
    --polyphen=b \\
    --custom "${gnomad},gnomADg,vcf,exact,0,AC,AC_afr,AC_afr_female,AC_afr_male,AC_ami,AC_ami_female,AC_ami_male,AC_amr,AC_amr_female,AC_amr_male,AC_asj,AC_asj_female,AC_asj_male,AC_eas,AC_eas_female,AC_eas_male,AC_female,AC_fin,AC_fin_female,AC_fin_male,AC_male,AC_nfe,AC_nfe_female,AC_nfe_male,AC_oth,AC_oth_female,AC_oth_male,AC_raw,AC_sas,AC_sas_female,AC_sas_male,AF,AF_afr,AF_afr_female,AF_afr_male,AF_ami,AF_ami_female,AF_ami_male,AF_amr,AF_amr_female,AF_amr_male,AF_asj,AF_asj_female,AF_asj_male,AF_eas,AF_eas_female,AF_eas_male,AF_female,AF_fin,AF_fin_female,AF_fin_male,AF_male,AF_nfe,AF_nfe_female,AF_nfe_male,AF_oth,AF_oth_female,AF_oth_male,AF_raw,AF_sas,AF_sas_female,AF_sas_male,nhomalt,nhomalt_afr,nhomalt_afr_female,nhomalt_afr_male,nhomalt_ami,nhomalt_ami_female,nhomalt_ami_male,nhomalt_amr,nhomalt_amr_female,nhomalt_amr_male,nhomalt_asj,nhomalt_asj_female,nhomalt_asj_male,nhomalt_eas,nhomalt_eas_female,nhomalt_eas_male,nhomalt_female,nhomalt_fin,nhomalt_fin_female,nhomalt_fin_male,nhomalt_male,nhomalt_nfe,nhomalt_nfe_female,nhomalt_nfe_male,nhomalt_oth,nhomalt_oth_female,nhomalt_oth_male,nhomalt_raw,nhomalt_sas,nhomalt_sas_female,nhomalt_sas_male" \\
    -o TMP/jeter.vcf \\
    --quiet \\
    --vcf \\
    --no_stats \\
    --plugin SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel} \\
    --plugin LOEUF,file=${loeuf},match_by=gene \\
    --plugin UTRAnnotator,file=${utrfile} 

bcftools view --threads ${task.cpus} -O z -o "${prefix}.vcf.gz" TMP/jeter.vcf
bcftools index --threads ${task.cpus} -t -f "${prefix}.vcf.gz"

cat << END_VERSIONS > versions.yml
"${task.process}":
    db: todo
END_VERSIONS
"""
}
