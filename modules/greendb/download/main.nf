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
process GREENDB_DOWNLOAD {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.bed.gz"),path("*.bed.gz.tbi"),path("*.header"),emit:bed
    path("versions.yml"),emit:versions
script:
    def TAG = task.ext.tag?:"GREENDB"
    def url = task.ext.url?:""
    if(url.isEmpty() && meta.ucsc_name!=null) {
        if(meta.ucsc_name == "hg38") {
            url = "https://zenodo.org/record/5636209/files/GRCh38_GREEN-DB.bed.gz?download=1";
        } else if(meta.ucsc_name == "hg19") {
            url = "https://zenodo.org/record/5636209/files/GRCh37_GREEN-DB.bed.gz?download=1";
        } 
    }

    def whatis = "GREEN-DB is a comprehensive collection of 2.4 million regulatory elements in the human genome collected from previously published databases, high-throughput screenings and functional studies. ${url}"
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP"
    def prefix = task.ext.prefix?:"${meta.id}.${TAG}"
"""
hostname 1>&2
mkdir -p TMP
env | grep -i proxy 1>&2

if ${!url.isEmpty()}
then

    curl -o TMP/jeter.bed.gz -L "${url}"

    gunzip -c TMP/jeter.bed.gz |\\
        grep -v '^chromosome' |\\
        cut -f1-3,5 |\\
        jvarkit ${jvm} bedrenamechr -f "${dict}" --column 1 --convert SKIP  |\\
        LC_ALL=C sort  --buffer-size=${task.memory.mega}M  -T . -t '\t' -k1,1 -k2,2n |\\
        uniq |\\
        bgzip > TMP/${prefix}.bed.gz


    echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis}. ${url}.">' > ${prefix}.header


else
    echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="GREENDB is not available">' > ${prefix}.header
    touch  TMP/${prefix}.bed
    bgzip  TMP/${prefix}.bed
fi


tabix --force -p bed TMP/${prefix}.bed.gz

mv TMP/${prefix}.bed.gz ./
mv TMP/${prefix}.bed.gz.tbi ./



cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
 def prefix = task.ext.prefix?:"${meta.id}.GREENDB"
"""
touch ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.header versions.yml
"""
}
