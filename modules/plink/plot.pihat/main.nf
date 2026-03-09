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
process PLOT_PIHAT {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta), path(genome)
output:
    tuple val(meta),path("*.png"),optional:true,emit:png
    tuple val(meta),path("*.pdf"),optional:true,emit:pdf
    tuple val(meta),path("*_mqc.txt"),optional:true,emit:high_txt
    path("versions.yml"),emit:versions
script:
    def format = task.ext.format?:(meta.format?:"png")
    def prefix=task.ext.prefix?:"${meta.id}.pihat"
    def max_pihat = task.ext.max_pihat?:0.1
    def plot_size = task.ext.plot_size?:1000

"""
mkdir -p TMP

# check input file
awk '(NR==1) {print \$10;}' '${genome}' | grep PI_HAT

cat << '__EOF__' > TMP/jeter.R
genome <- read.table(file="${genome}",header=TRUE,stringsAsFactors=FALSE)

head(genome)

${format}("TMP/${prefix}_mqc.${format}",width = ${plot_size}, height = ${plot_size}, unit = "px")
plot(genome\$PI_HAT,
    ylim=c(0,1.0),
    xlab="Individuals Pair",
    ylab="PI-HAT",
    main="PI-HAT",
    sub="pihat for each pair of sample"
)
abline(h=${max_pihat},col="blue");
dev.off()
__EOF__


R --no-save < TMP/jeter.R

mv -v TMP/${prefix}* ./

awk 'NR==1 || ((\$10 * 1.0) > ${max_pihat})' '${genome}' > '${prefix}.high_pihat_mqc.txt'


cat << EOF > versions.yml
${task.process}:
    R: \$(R --version | awk '(NR==1) {print \$3;}')
EOF
"""

stub:
    def format = task.ext.format?:(meta.format?:"png")
    def prefix=task.ext.prefix?:"${meta.id}.pihat"
"""
touch versions.yml '${prefix}.${format}' '${prefix}.high_pihat_mqc.txt'
"""
}