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

process ENSEMBL_REG_DOWNLOAD{
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
   	def TAG = task.ext.tag?:"ENSEMBL_REG"
	def whatis="Ensembl regulatory Features"
    def url= task.ext.url?:"";
    if(url.isEmpty() && meta.ucsc_name!=null) {
            if(meta.ucsc_name == "hg38") {
                url = "https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v114.gff3.gz"
            } else if(meta.ucsc_name == "hg19") {
                url = "https://ftp.ensembl.org/pub/release-114/regulation/homo_sapiens/GRCh37/annotation/Homo_sapiens.GRCh37.regulatory_features.v114.gff3.gz"
            } 
        }
    def prefix = task.ext.prefix?:"${meta.id}.${TAG}"
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
mkdir -p TMP

if ${!url.isEmpty()}
then

    curl -L "${url}" |\\
            gunzip -c |\\
            awk -F '\t' '/^[^#]/{printf("%s\t%d\t%s\t%s\\n",\$1,int(\$4)-1,\$5,\$3);}' |\\
            jvarkit ${jvm} bedrenamechr -R ${dict} --column 1 --convert SKIP  |\\
            LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k1,1 -k2,2n -T TMP |\\
            uniq | bgzip > TMP/${prefix}.bed.gz
        

        echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="${whatis} ${url}">' > ${prefix}.header

else


    echo '##INFO=<ID=${TAG},Number=.,Type=String,Description="ENSEMBL REG not available">' > ${prefix}.header

    touch TMP/${prefix}.bed
    bgzip TMP/${prefix}.bed

fi

tabix -p bed -f TMP/${prefix}.bed.gz

mv TMP/${prefix}.bed.gz ./
mv TMP/${prefix}.bed.gz.tbi ./


cat << EOF > versions.yml
${task.process}:
	url: "${url}"
EOF
"""


stub:
    def prefix = "${meta.id}.ENSREG"
"""
touch ${prefix}.bed.gz ${prefix}.bed.gz.tbi ${prefix}.header versions.yml
"""
}