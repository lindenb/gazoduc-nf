/*

Copyright (c) 2026 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


baserefseq="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene"

include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


workflow BUILD_NCBI_GENE_DB {
    take:
        meta
    main:
	version_ch = Channel.empty()

        awk_ch = MAKE_AWK_SCRIPT(meta)
	version_ch = version_ch.mix(awk_ch.version)
    
        installed_ch = WGET_INSTALLED(meta)
	version_ch = version_ch.mix(installed_ch.version)


        gbff_ch = installed_ch.output.splitCsv(sep:'\t',header:false).
            map{T->T[1]}.
            filter{T->T.endsWith(".gbff.gz")}

        gbff2tsv_ch = WGET_GBFF(meta, awk_ch.output, gbff_ch)
	version_ch = version_ch.mix(gbff2tsv_ch.version)

        merge_ch = MERGE_GBFF(meta, gbff2tsv_ch.output.collect())
	version_ch = version_ch.mix(merge_ch.version)

	version_ch = MERGE_VERSION(meta, "ncbigene", "NCBI genes", version_ch.collect())

    emit:
	version= version_ch
	output = merge_ch.output
    }

process WGET_INSTALLED {
executor "local"
input:
    val(meta)
output:
    path("refseqgene.files.installed"),emit:output
    path("version.xml"),emit:version
script:
"""
wget -O refseqgene.files.installed "${baserefseq}/refseqgene.files.installed"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download installed</entry>
</properties>
EOF
"""
}


process MAKE_AWK_SCRIPT {
executor "local"
input:
    val(meta)
output:
    path("gb2tsv.awk"),emit:output
    path("version.xml"),emit:version
"""

cat << '__EOF__' > gb2tsv.awk
/^VERSION / {LOCUS=substr(\$0,13);next;}
/^COMMENT / {state=1;DESC=substr(\$0,13);next;}
/^\\/\\// {printf("%s\t%s\\n",LOCUS,DESC);LOCUS=""; DESC="";state=0;next;}
/^[^ ]/ {state=0;next;}
/^[ ]/ {if(state==1) DESC=sprintf("%s%s",DESC,substr(\$0,13));}   
__EOF__


#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">awk script parsing genbank</entry>
</properties>
EOF
"""
}

process WGET_GBFF {
tag "${gbff}"
maxForks 6
afterScript "rm -rf TMP"
input:
    val(meta)
    path(awk)
    val(gbff)
output:
    path("gbff.tsv.gz"),emit:output
    path("version.xml"),emit:version
script:
	def url = "${baserefseq}/${gbff}"
"""
set -o pipefail

mkdir -p TMP

wget -O - "${url}" |\
    gunzip -c |\
    awk -f "${awk}" |\
    sed -e 's/\\t.*Summary:[ ]*/\\t/' -e 's/\\.[ ]*\\[provided[ ]*by[ ]*RefSeq,[ A-Za-z0-9]*\\]\\./\\./i' -e 's/REVIEWED REFSEQ: This record has been curated by NCBI staff\\.//i'  |\
    awk -F '\t' '{OFS="\t";gsub(/\\.[0-9]*\$/,"",\$1); print;}' |\
    LC_ALL=C sort -T TMP -t '\t' -k1,1 > TMP/gbff.tsv

gzip --best TMP/gbff.tsv

mv TMP/gbff.tsv.gz ./

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download gbff ${gbff}</entry>
	<entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

process MERGE_GBFF {
tag "N=${L.size()}"
//afterScript "rm -rf TMP"
input:
    val(meta)
    val(L)
output:
	path("${meta.prefix?:""}ncbi.genes.refseq.tsv.gz"),emit:output
	path("version.xml"),emit:version
script:
"""
set -o pipefail
mkdir -p TMP

wget -O - "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/Homo_sapiens.gene_info.gz" | gunzip -c | cut -f 2,3,9,10 | LC_ALL=C  sort -T TMP -t '\t' -k1,1 | uniq > TMP/jeter.a

# sort on geneid, remote version from accession number
wget -O - "${baserefseq}/LRG_RefSeqGene" | cut -f 2,3,4 | LC_ALL=C  sort -T TMP -t '\t' -k1,1 | awk -F '\t' '{OFS="\t";gsub(/\\.[0-9]*\$/,"",\$3); print;}' | uniq > TMP/jeter.b

LC_ALL=C join -t '\t' -1 1 -2 1 -o '1.2,1.1,1.3,1.4,2.3'  -a 1 -a 2 -e '.' TMP/jeter.a TMP/jeter.b  | LC_ALL=C  sort -T TMP -t '\t' -k5,5 > TMP/jeter.c

# merge, remove version from accession number
gunzip -c ${L.join(" ")} |  awk -F '\t' '{OFS="\t";gsub(/\\.[0-9]*\$/,"",\$1); print;}'  | LC_ALL=C sort -T TMP -t '\t' -k1,1 --unique > TMP/jeter.d

echo "gene_name\tgene_id\tgene_title\tgene_biotype\trefSeq\tgene_description" > TMP/nci.genes.refseq.tsv

LC_ALL=C join -t '\t' -1 5 -2 1 -o '1.1,1.2,1.3,1.4,1.5,2.2'  -a 1 -a 2 -e '.' TMP/jeter.c TMP/jeter.d |\
    LC_ALL=C  sort -T TMP -t '\t' -k2,2 --unique |\
    LC_ALL=C  sort -T TMP -t '\t' -k1,1 >> TMP/ncbi.genes.refseq.tsv

gzip --best TMP/ncbi.genes.refseq.tsv

mv TMP/ncbi.genes.refseq.tsv.gz "./${meta.prefix?:""}ncbi.genes.refseq.tsv.gz"


#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Join gbff</entry>
</properties>
EOF

"""
}
