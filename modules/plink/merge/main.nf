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
process PLINK_MERGE_BIM_BED_FAM {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"

input:
    tuple val(meta),path("PLINK/*")
output:
	tuple val(meta),path("*.bim"),path("*.bed"),path("*.fam"),emit:plink
    path("versions.yml")
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def args1  = "--const-fid 1 --allow-extra-chr --allow-no-sex "

"""
mkdir -p TMP

find PLINK/  -name "*.bim" | sed 's/\\.bim\$//' | sort -V -T TMP > TMP/jeter.list
test -s TMP/jeter.list

plink \\
    ${args1} \\
    --threads ${task.cpus} \\
    --merge-list TMP/jeter.list \\
    --make-bed \\
    --out ${prefix}

touch versions.yml
"""
stub:
"""
touch versions.yml ${prefix}.bim ${prefix}.bed ${prefix}.fam
"""
}
