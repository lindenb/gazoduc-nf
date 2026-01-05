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
include {verify     } from '../../../modules/utils/functions.nf'
include {isBlank    } from '../../../modules/utils/functions.nf'

process FASTQ_MERGE {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(fastqs)
output:
    tuple val(meta),path("*.fastq.gz"),emit:fastq
    path("versions.yml"),emit:versions
script:
    def L0 = fastqs
    if(fastqs instanceof Path) {
        log.warn("FASTQ_MERGE used with only one fastq ? ${fastqs}")
        L0 = [fastqs]
        }
  
    def prefix  = task.ext.prefix?:"${meta.prefix}" //yes,prefix, not ID
    verify(!isBlank(prefix),"${task.process} prefix is blank")
"""
mkdir -p TMP

cat ${L0.collect{it.name}.sort().join(" ")} > TMP/${prefix}.fastq.gz
mv TMP/${prefix}.fastq.gz ./

cat << EOF > versions.yml
${task.process}:
    cat: "\$( cat --version | head -n1 | awk '{print \$NF}' )"
EOF
"""

stub:
    def prefix  = task.ext.prefix?:"${meta.prefix}" //yes,prefix, not ID
    verify(!isBlank(prefix),"${task.process} prefix is blank")
"""
touch versions.yml ${prefix}.fastq.gz
"""
}