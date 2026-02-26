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
afterScript "rm -rf TMP"
label "process_short" //single is not enough memory
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
    tuple val(meta),path(bim),path(bed),path(fam)
output:
    tuple val(meta),path("*.genome"),emit:genome
    tuple val(meta),path("merged.*"),emit:merged_plink
    path("versions.yml"),emit:versions
script:
    def num_mds_components = task.ext.mds_components?:10
    def plink_args  = "--const-fid 1 --allow-extra-chr --allow-no-sex --threads ${task.cpus}"
    def prefix = task.ext.prefix?:meta.id?:"pihat"
"""
mkdir -p TMP
set -x

## IBS matrix
plink --bfile ${bim.baseName} --genome --out TMP/matIBS_data
wc -l TMP/matIBS_data.genome 1>&2
cp TMP/matIBS_data.genome "${prefix}.genome"


## Remove related individuals
awk '(NR==1 || \$10>0.1) {printf("%s%s\t%s\\n",(NR==1?"#":""),\$3 ,\$4);}' TMP/matIBS_data.genome |\\
    LC_ALL=C sort |\\
    uniq  |\\
    sed 's/^#//' > TMP/indiv_pihat_sup0.1.txt
wc -l TMP/indiv_pihat_sup0.1.txt 1>&2

cp TMP/indiv_pihat_sup0.1.txt ./

plink \\
    ${plink_args} \\
    --bfile TMP/merged \\
    --remove TMP/indiv_pihat_sup0.1.txt \\
    --make-bed \\
    --out TMP/data_sansApp

## MDS
plink \\
    ${plink_args} \\
    --bfile TMP/data_sansApp \\
    --read-genome TMP/matIBS_data.genome \\
    --mds-plot ${num_mds_components} \\
    --cluster \\
    --out TMP/strat_data

mv TMP/strat_data.mds ${prefix}.mds


mv TMP/merged.* ./

cat << EOF > versions.yml
${task.process}:
    plink: "\$(plink --version | awk '{print \$2}')"
EOF
"""
}



