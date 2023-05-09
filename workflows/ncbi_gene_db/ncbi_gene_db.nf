params.baserefseq="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene"
workflow {
    BUILD_NCBI_GENE_DB(params)
    }

workflow BUILD_NCBI_GENE_DB {
    take:
        meta
    main:
        awk_ch = MAKE_AWK_SCRIPT(meta)

        gbff_ch = WGET_INSTALLED(meta).output.splitCsv(sep:'\t',header:false).
            map{T->T[1]}.
            filter{T->T.endsWith(".gbff.gz")}

        gbff2tsv_ch = WGET_GBFF(meta, awk_ch.output, gbff_ch)

        MERGE_GBFF(meta, gbff2tsv_ch.output.collect())

    }

process WGET_INSTALLED {
executor "local"
input:
    val(meta)
output:
    path("refseqgene.files.installed"),emit:output
script:
"""
wget -O refseqgene.files.installed "${meta.baserefseq}/refseqgene.files.installed"
"""
}


process MAKE_AWK_SCRIPT {
executor "local"
input:
    val(meta)
output:
    path("gb2tsv.awk"),emit:output
"""

cat << '__EOF__' > gb2tsv.awk
/^VERSION / {LOCUS=substr(\$0,13);next;}
/^COMMENT / {state=1;DESC=substr(\$0,13);next;}
/^\\/\\// {printf("%s\t%s\\n",LOCUS,DESC);LOCUS=""; DESC="";state=0;next;}
/^[^ ]/ {state=0;next;}
/^[ ]/ {if(state==1) DESC=sprintf("%s%s",DESC,substr(\$0,13));}   
__EOF__

"""
}

process WGET_GBFF {
executor "local"
tag "${gbff}"
maxForks 1
input:
    val(meta)
    path(awk)
    val(gbff)
output:
    path("gbff.tsv"),emit:output
script:
"""
set -o pipefail

wget -O - "${meta.baserefseq}/${gbff}" |\
    gunzip -c |\
    awk -f "${awk}" |\
    sed 's/\\t.*Summary:[ ]*/\\t/' |\
    sed 's/\\.[ ]*\\[provided by[ ]*RefSeq,[ A-Za-z0-9]*\\]\\./\\./' |\
    LC_ALL=C sort -T . -t '\t' -k1,1 > gbff.tsv
"""
}

process MERGE_GBFF {
tag "N=${L.size()}"
afterScript "rm -f jeter.a jeter.b"
input:
    val(meta)
    val(L)
script:
"""
set -o pipefail

wget -O - "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/gene_RefSeqGene" | cut -f 2- | LC_ALL=C  sort -T . -t '\t' -k3,3 > jeter.a

LC_ALL=C sort -T . -t '\t' -k1,1 --merge ${L.join(" ")} > jeter.b

LC_ALL=C join -t '\t' -1 3 -2 1 -o '1.2,1.1,1.3,2.2'  -a 1 -a 2 -e '.' jeter.a jeter.b |\
    LC_ALL=C  sort -T . -t '\t' -k1,1 |\
    gzip --best > 20230509.nci.genes.refseq.tsv.gz

"""
}
