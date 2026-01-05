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
include {isBlank      } from  '../../../modules/utils/functions.nf'

/** download blacklisted regions from encode */
process ENCODE_BLACKLIST {
tag "${meta1.id?:fasta.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(dict)
output:
	tuple val(meta1),path("*.bed"),emit:bed
	path("versions.yml"),emit:versions
script:
	def url = task.ext.url?:"";
    def prefix = task.ext.suffix?:""
    if(!isBlank(url)) {
        //nothing
        }
    else if(isBlank(meta1.ucsc_name)) {
        url = ""
        }
    else if(meta1.ucsc_name == "hg19") {
        url = "https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz";
        if(isBlank(prefix)) prefix = "encode.hg19.blacklist.ENCFF001TDO"
        }
    else if(meta1.ucsc_name == "hg38") {
        url = "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz";
        if(isBlank(prefix)) prefix = "encode.hg38.blacklist.ENCFF356LFX"
        }
    else {
        url = ""
        }
    if(isBlank(prefix)) prefix="encode.blacklist";
    def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
    if(isBlank(url)) {
        log.warn("Empty url for ${task.process}");
        }
"""
hostname 1>&2
mkdir -p TMP

if ${!isBlank(url)}
then
    curl -L -o TMP/jeter.out  "${url}"
    if  (file TMP/jeter.out | grep "gzip compressed")
    then
        gunzip -c  TMP/jeter.out > TMP/jeter.bed
    else
        mv  TMP/jeter.out TMP/jeter.bed
    fi

else

    touch TMP/jeter.bed

fi

cut -f1,2,3 TMP/jeter.bed |\\
    awk '/^#/ {next;} {print;}' |\\
	jvarkit  bedrenamechr --column 1 -R "${dict}" --convert SKIP |\\
    LC_ALL=C sort --buffer-size=${task.memory.mega}M -t '\t' -k1,1 -k2,2n -T TMP |\\
    bedtools merge > TMP/jeter2.bed

mv TMP/jeter2.bed "${prefix}.bed"

cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
END_VERSIONS
"""

stub:
"""
touch versions.yml blacklist.bed
"""
}
