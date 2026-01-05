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
process DOWNLOAD_UTR_ANNOTATOR {
tag "${meta.id}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(dict)
output:
    tuple val(meta),path("*.txt"),emit:output
    path("versions.yml"),emit:versions
script:
    def filename= task.ext.filename?:"uORF_5UTR_undefined_PUBLIC.txt"
    def url = task.ext.url?:""
    if(url.isEmpty()) {
        if((meta.ucsc_name==null)) {
            throw new IllegalArgumentException("${task.process} missing ucsc_name");
        } else if(meta.ucsc_name =="hg38") {
            filename="uORF_5UTR_GRCh38_PUBLIC.txt"
        } else if(meta.ucsc_name == "hg19") {
            filename="uORF_5UTR_GRCh37_PUBLIC.txt"
        } else {
            throw new IllegalArgumentException("${task.process} unknown ucsc_name");
        }
        url = "https://github.com/ImperialCardioGenetics/UTRannotator/raw/refs/heads/master/${filename}"
        }
"""
curl -L -o "${filename}" "${url}"

cat << END_VERSIONS > versions.yml
"${task.process}":
    url: "${url}"
END_VERSIONS
"""

stub:
"""
touch uORF_5UTR_undefined_PUBLIC.txt versions.yml
"""
}

