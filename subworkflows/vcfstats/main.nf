include {BCFTOOLS_STATS} from '../../modules/bcftools/stats'
workflow VCF_STATS {
take:
    meta
    fasta
    fai
    dict
    gtf
    gff3
    pedigree
    sample2pop
    bed
    vcfs
main:
    versions= Channel.empty()
    multiqc = Channel.empty()


    BCFTOOLS_STATS(
        fasta,
        fai,
        bed,
        gtf,
        [[id:"no_samples"],[]],
        vcfs
        )
    versions= versions.mix(BCFTOOLS_STATS.out.versions)

emit:
    versions
    multiqc
}

process F_MISSING {
input:
    tuple val(meta),path("VCFS")
script:
"""
mkdir -p TMP
find VCFS/ \\( -name "*.bcf" -o -name "*.vcf.gz" \\) | sort > TMP/jeter.list
bcftools concat --file-list TMP/jeter.list
"""
}