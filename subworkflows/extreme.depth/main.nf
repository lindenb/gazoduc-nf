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

include {DEPTH_OUTLIERS         } from '../../modules/samtools/depth.outliers'

workflow EXTREME_DEPTH {
 take:
    metadata
    fasta
    fai
    dict
    bed
    bams // [meta,bam,bai]
main:
    versions = Channel.empty()
    multiqc = Channel.empty()
   
    DEPTH_OUTLIERS( fasta, fai, bed, bams  )
    versions = versions.mix(DEPTH_OUTLIERS.out.versions)

    
    
emit:
    versions
    multiqc
}


process MERGE_OUTLIERS {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(bed)
    tuple val(meta ),path("BED/")
output:
    tuple val(meta ),path("*.outliers.bed"),emit:outliers
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.outliers"
    def args1=task.ext.args1?:"-d 2"
"""
mkdir -p TMP
find BED/ -name "*.bed.gz"  | while read B
do
    gunzip -c "\${B}" | cut -f1,2,3 | sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter1.bed
    if test ! -f TMP/jeter2.bed
    then
       mv TMP/jeter1.bed TMP/jeter2.bed
    else
        bedtools intersect -sorted -a TMP/jeter1.bed -b TMP/jeter2.bed |\\
            sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter3.bed
        mv TMP/jeter3.bed  TMP/jeter2.bed
    fi
done

bedtools merge ${args1} -i TMP/jeter2.bed > TMP/jeter.bed

if  test -s TMP/jeter.bed
then
    bedtools intersect -a bed -b TMP/jeter.bed |\\
        sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter2.bed
    mv TMP/jeter2.bed ${prefix}.adjusted.bed




mv TMP/jeter.bed ${prefix}.outliers.bed

"""
}