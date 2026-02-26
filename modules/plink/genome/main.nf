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
process PLINK_GENOME {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta),path(bim),path(bed),path(fam)
output:
	tuple val(meta),path("*.genome"),emit:genome
    tuple val(meta),path("*.related.tsv"),emit:related
    tuple val(meta),path("*.log"),emit:log
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def args1  = "--const-fid 1 --allow-extra-chr --allow-no-sex "
    def treshold = task.ext.treshold?:"0.1"
"""
mkdir -p TMP

plink \\
    --threads ${task.cpus} \\
    --bfile ${bim.baseName} \\
    --genome \\
    --out TMP/${prefix}


## Extract related individuals
awk '(NR==1 || \$10 > ${treshold}) {printf("%s%s\t%s\\n",(NR==1?"#":""),\$3 ,\$4);}' TMP/${prefix}.genome |\\
    LC_ALL=C sort -T TMP |\\
    uniq  |\\
    sed 's/^#//' > TMP/${prefix}.related.tsv

mv -v TMP/${prefix}* ./

cat << EOF > versions.yml
${task.process}:
    plink: \$(plink --version)
EOF
"""
stub:
"""
touch versions.yml ${prefix}.genome
"""
}
