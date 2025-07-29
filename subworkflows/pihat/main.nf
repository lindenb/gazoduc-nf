include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed/main.nf'

workflow PIHAT {
    take:
        meta
        fasta
        fai
        dict
        vcfs
    main:


        versions = Channel.empty()
        
        VCF_TO_CONTIGS(vcfs)
        versions = versions.mix(VCF_TO_CONTIGS.out.versions)

        VCF_TO_CONTIGS.out.output
            .splitCsv(header:false,sep:'\t',elem:1)
            .map{it[0]}
            .filter(it.matches("(chr)?[0-9]+"))
            .map{[[id:it],it]}
        


        if(fasta[0].ucsc_name && fasta[0].ucsc_name.equals("hg38")) {
            ch1 = each_contig.map{[
                it[0],
                it[1],
                "/LAB-DATA/GLiCID/projects/BiRD_resources/species/human/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.${it[1].startsWith("chr")?it[1]:"chr"+it[1]}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
                ]}
            .map{[
                it[0],
                it[1],
                file(it[2]),
                file(it[2]+".tbi")
                ]}
            .filter{it[2].exists() && it[3].exists()}
            }
        else
            {
            ch1 = Channel.empty()
            }

    emit:
        versions
}