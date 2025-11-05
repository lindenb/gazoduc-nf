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
process DRAGMAP_INDEX {
tag "${fasta.name}"
label "process_medium"
conda "${moduleDir}/../../../conda/dragmap.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta ),path(fasta)
    tuple val(meta2),path(bed_mask)
output:
	tuple val(meta),path("DragenIndex_*"),emit:dragmap_index
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args = task.ext.args?:""
"""
    mkdir -p TMP/${prefix}

    dragen-os \\
        ${args} \\
        --ht-num-thread ${task.cpus} \\
        --build-hash-table true \\
        --ht-reference ${fasta}  \\
        --output-directory TMP/${prefix} \\
        --output-file-prefix=${prefix} \\
        ${bed_mask?"--ht-mask-bed=${bed_mask}":""}
    

    mv -v TMP/${prefix} "DragenIndex_${prefix}"

cat << END_VERSIONS > versions.yml
${task.process}:
    dragmap : \$(dragen-os -V)
END_VERSIONS
"""

stub:
"""
mkdir -p TMP
(cd TMP && touch hash_table.cfg hash_table.cfg.bin hash_table.cmp hash_table_stats.txt reference.bin ref_index.bin repeat_mask.bin str_table.bin )
mv TMP "DragenIndex_${fasta.baseName}"
touch versions.yml
"""
}
