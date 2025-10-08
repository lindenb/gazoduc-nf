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
include {k1_signature} from '../../utils/k1.nf'

process BCFTOOLS_GUESS_PLOIDY {
label "process_single"
tag "${vcf.name} ${meta.id?:""}"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
    tuple val(meta1),path(fasta)
    tuple val(meta2),path(fai)
    tuple val(meta),path(vcf),path(idx)
output:
    tuple val(meta),path("*.ploidy.txt"), emit:output
    tuple val(meta),path("*.png"), optional:true, emit:image
    path("versions.yml"),emit:versions
when:
    task.ext.when == null || task.ext.when
script:
    def k1 = k1_signature()
    def prefix = task.ext.prefix?:(meta.id?:vcf.simpleName)
    def args = task.ext.args?:""
"""
mkdir -p TMP

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg38}\tb38
1:${k1.hg19}\tb37
chr1:${k1.hg38}\thg38
chr1:${k1.hg19}\thg19
EOF

awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
    sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
    sort |\\
    uniq > TMP/jeter.build

BUILD=`cat  TMP/jeter.build`
test ! -z "\${BUILD}"

bcftools +guess-ploidy ${args} ${vcf} -g "\${BUILD}" > ${prefix}.ploidy.txt

guess-ploidy.py "${prefix}.ploidy.txt" ${prefix}.ploidy.img || true

cat << END_VERSIONS > versions.yml
"${task.process}":
    bcftools: \$(bcftools version | awk '(NR==1)  {print \$NF}')
    build: \${BUILD}
END_VERSIONS
"""

stub:
"""
echo '${meta.id}\ttodo' > ${meta.id}.ploidy.txt
touch  ${meta.id}.ploidy.png
touch versions.yml
"""
}
