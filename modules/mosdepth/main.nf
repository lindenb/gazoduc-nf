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

include {parseBoolean} from '../../modules/utils/functions.nf'

process MOSDEPTH {
tag "${meta.id?:""}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/mosdepth.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta), path(bam), path(bai),path(optional_bed)
output:
        tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
        tuple val(meta), path('*.summary.txt')          , emit: summary_txt
        tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
        tuple val(meta), path('*.regions.bed.gz'), path('*.regions.bed.gz.csi'), optional:true, emit: regions_bed
        path  "versions.yml"                            , emit: versions
script:
    def per_base = parseBoolean(task.ext.per_base?:false)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP

mosdepth \\
    --threads ${task.cpus} \\
   ${per_base==true?"":"--no-per-base"} \\
   ${optional_bed?"--by \"${optional_bed}\"":""} \\
    --fasta "${fasta}" \\
    ${args} \\
    TMP/${prefix} \\
    ${bam}

	
mv -v TMP/${prefix}* ./

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
END_VERSIONS
"""

stub:
"""
touch versions.yml
touch ${meta.id}.global.dist.txt
echo "chrom\tmean" > ${meta.id}.summary.txt
touch ${meta.id}.region.dist.txt
touch ${meta.id}.region.bed.bed.gz
touch ${meta.id}.region.bed.bed.gz.csi
"""
}
