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
process RUN_ORAD {
label 'process_medium'
tag "${meta.id} ${ora_file.name}"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(oradir)
        tuple val(meta ),path(ora_file)
output:
        tuple val(meta),path("*q.gz",arity:'1..*'),emit:fastq
        path("versions.yml"),emit:versions
when:
        task.ext.when == null || task.ext.when
script:
"""
mkdir -p TMP

${oradir}/orad \
        --threads ${task.cpus} \\
        --path TMP \\
        --ora-reference "${oradir}/oradata" \\
        --name ${ora_file} \\
        -q

mv -v TMP/*.gz ./
    
cat <<-END_VERSIONS > versions.yml    
"${task.process}":
    orad: \$(${oradir}/orad --version | awk '(NR==1) {print \$3;}')
END_VERSIONS
"""

stub:
"""
touch ${meta.id}.R1.fq.gz
touch ${meta.id}.R2.fq.gz
touch versions.yml
"""
}
