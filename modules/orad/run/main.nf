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
include {removeCommonSuffixes  } from '../../utils/functions.nf'

process RUN_ORAD {
label 'process_medium'
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
        tuple val(meta1),path(oradir)
        tuple val(meta ),path(fastq_files)
output:
        tuple val(meta),path("*q.gz",arity:'1..*'),emit:fastq
        path("versions.yml"),emit:versions
script:
        def args1 =  task.ext.args1?:""
        def prefix = task.ext.prefix?:(fastq_files instanceof List && fastq_files.size()==1 ? removeCommonSuffixes(fastq_files[0].name):"${meta.id}")
"""
mkdir -p TMP



${oradir}/orad \\
        ${args1} \\
        --threads ${task.cpus} \\
        --path TMP \\
        --ora-reference "${oradir}/oradata" \\
        -q \\
        - < ${fastq_files}

if test -f TMP/-R1.fastq.gz
then
        mv TMP/-R1.fastq.gz TMP/${prefix}_R1.fastq.gz
fi

if test -f TMP/-R2.fastq.gz
then
        mv TMP/-R2.fastq.gz TMP/${prefix}_R2.fastq.gz
fi



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
