workflow PIHAT {
    take:
        meta
        fasta
        fai
        dict
        vcfs
    main:
        versions = Channel.empty()

    emit:
        versions
}