/*

Copyright (c) 2025 Pierre Lindenbaum

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
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(dict)
output:
    tuple val(meta),path("*.bed.gz"), path("*.bed.gz.tbi"), path("*.header"),emit:bed
    path("versions.yml"),emit:versions
script:
    def base = task.ext.base?:"https://remap.univ-amu.fr/storage/remap2022/"
    def url= task.ext.url?:""
    if(url.isEmpty() && meta.ucsc_name!=null) {
        if(meta.ucsc_name == "hg19") {
            url ="${base}/hg19/MACS2/remap2022_crm_macs2_hg19_v1_0.bed.gz"
        } else  if(meta.ucsc_name == "hg38") {
            url ="${base}/hg38/MACS2/remap2022_crm_macs2_hg38_v1_0.bed.gz"
        }
    }
    
    def TAG = task.ext.tag?:"REMAP"
    def WHATIZ = "ReMap is a large scale integrative analysis of DNA-binding experiments for Homo sapiens, Mus musculus, Drosophila melanogaster and Arabidopsis thaliana transcriptional regulators. The catalogues are the results of the manual curation of ChIP-seq, ChIP-exo, DAP-seq from public sources (GEO, ENCODE, ENA)."
    def prefix = task.ext.prefix?:"${meta.id}.${TAG}"
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

if ${!url.isEmpty()}
then

curl -L "${url}"  |\\
    gunzip -c |\\
    cut -f1,2,3 |\\
    jvarkit ${jvm} bedrenamechr -R "${dict}" --column 1 --convert SKIP  |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k1,1 -k2,2n -T TMP |\\
    bedtools merge |\\
    sed 's/\$/\t1/' |\\
    bgzip > TMP/${prefix}.bed.gz

    echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="${WHATIZ}">' > ${prefix}.header

else

    echo '##INFO=<ID=${TAG},Number=0,Type=Flag,Description="Remap is not available for this organism/build.">' > ${prefix}.header

    touch TMP/${prefix}.bed
    bgzip TMP/${prefix}.bed
fi

tabix -p bed -f TMP/${prefix}.bed.gz

mv TMP/${prefix}.bed.gz ./
mv TMP/${prefix}.bed.gz.tbi ./



cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}.REMAP"

"""
touch {prefix}.bed.gz {prefix}.bed.gz.tbi {prefix}.header versions.yml
"""
}
