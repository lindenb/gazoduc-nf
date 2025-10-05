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
process FASTP {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/fastp.yml"
input:
    tuple val(meta),path(fastqs)
output:
    tuple val(meta), path('*.fq.gz', arity: '1..*') , emit: fastqs
    tuple val(meta), path('*.json') ,emit: json
    tuple val(meta), path('*.html') ,emit: html
    path("versions.yml"),emit:versions
script:
    if(!(fastqs instanceof List)) throw new IllegalArgumentException("${task.process}: fastqs should be a List");
    if(fastqs.size()>2) throw new IllegalArgumentException("${task.process}: fastqs.size > 2");
    def L = fastqs.sort()
    def prefix = task.ext.prefix?:""
    def args1 = task.ext.args1?:""
"""
    mkdir -p TMP

    if ${L.size()==2}
    then
        fastp \\
             ${args1} \\
            --thread ${task.cpus} \\
            --json TMP/jeter.json \\
            --html TMP/jeter.html \\
            --report_title "${meta.id}" \\
            -i "${fastqs[0]}" \\
            -I "${fastqs[1]}" \\
            -o "TMP/${prefix}${L[0].baseName}.fastp.fq.gz" \\
            -O "TMP/${prefix}${L[1].baseName}.fastp.fq.gz"
    else
         fastp \\
            ${args1} \\
            --thread ${task.cpus} \\
            --json TMP/jeter.json \\
            --html TMP/jeter.html \\
            --report_title "${meta.id}" \\
            -i "TMP/${prefix}${L[0]}" \\
            -o "TMP/${prefix}${L[0].baseName}.fastp.fq.gz"
    fi

mv -v TMP/*.gz ./
mv TMP/jeter.json "${meta.id}.json"
mv TMP/jeter.html "${meta.id}.html"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
END_VERSIONS
"""
}