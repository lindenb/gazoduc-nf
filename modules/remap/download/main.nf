/*

Copyright (c) 2024 Pierre Lindenbaum

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
process REMAP_DOWNLOAD {
tag "${meta.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
when:
    task.ext.when == null || task.ext.when
input:
    tuple val(meta ),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
    tuple val(meta),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:bed
    path("versions.yml"),emit:versions
script:
    def base = "https://remap.univ-amu.fr/storage/remap2022/"
    if(!meta.ucsc_name || meta.ucsc_name.isEmpty()) throw new IllegalArgumentException("${task.process} ucsc mssing");
    def url="";
    if(meta.ucsc_name.equals("hg19")) {
        url ="${base}/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz"
    } else  if(meta.ucsc_name.equals("hg38")) {
        url ="${base}/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"
    } else {
         throw new IllegalArgumentException("${task.process} bad ucsc ");
    }
    
    def TAG = "REMAP"
    def WHATIZ = "ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA)."
"""
hostname 1>&2
mkdir -p TMP

curl -L "${url}"  |\\
    gunzip -c |\\
    cut -f1,2,3 |\\
    jvarkit -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP bedrenamechr -R "${fasta}" --column 1 --convert SKIP  |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k1,1 -k2,2n -T TMP |\\
    bedtools merge |\\
    sed 's/\$/\t1/' |\\
    bgzip > TMP/${TAG}.bed.gz


tabix -p bed -f TMP/${TAG}.bed.gz

mv TMP/${TAG}.bed.gz ./
mv TMP/${TAG}.bed.gz.tbi ./


echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="${WHATIZ}">' > ${TAG}.header

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""
}
