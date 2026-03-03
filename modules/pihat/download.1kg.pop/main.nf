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
process DOWNLOAD_1KG_SAMPLE2POP {
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
label "process_single"
input:
    val(meta)
output:
    tuple val(meta),path("*.tsv"),emit:tsv
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"integrated_call_samples_1kg"
    def url = task.ext.url?:"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20250704.ALL.ped"
"""
curl -L "${url}"  |\\
	sed '1s/[ ]/_/g' > "${prefix}.tsv"

cat << EOF > versions.yml
${task.process}:
    url: "${url}"
EOF
"""

stub:
 def prefix = task.ext.prefix?:"integrated_call_samples_1kg"
"""
touch ${prefix}.tsv versions.yml
"""
}
