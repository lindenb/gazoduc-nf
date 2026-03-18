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
 
process QQMAN {
label "process_single"
afterScript "rm -rf TMP"
tag "${meta.id}"
conda "${moduleDir}/../../conda/qqman.yml"
input:
	tuple val(meta ),path(table)
output:
    tuple val(meta),path("*.manhattan_mqc.*"),optional:true,emit:manhattan
    tuple val(meta),path("*.qqplot_mqc.*"),optional:true,emit:qqplot
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def format = task.ext.format?:"png"
    def title = task.ext.title?:"${meta.id}"
    def subtitle = task.ext.subtitle?:""
    def plot_size = task.ext.plot_size?:500
    def args1= task.ext.args1?:""
    def args2= task.ext.args2?:""
    def cmd = task.ext.cmd?:"cat" //give a chance to transforme the table
"""
mkdir -p TMP

if ${!cmd.isEmpty()}
then
	${table.name.endsWith(".gz")?"gunzip -c":"cat"} ${table} | ${cmd} > TMP/jeter.tsv
fi

cat << '__EOF__' > TMP/jeter.R
library("qqman")
T1 <- read.table("${cmd.isEmpty()?"${table}":"TMP/jeter.tsv"}",header=TRUE,sep="\t",stringsAsFactors=FALSE)
head(T1)
T1 <- T1[grepl("(chr)?[0-9XY]*", T1\$CHR), ]
head(T1)
T1\$CHR <- sub("chr", "", T1\$CHR)
head(T1)
T1\$CHR[T1\$CHR == "X"] <- 23
head(T1)
T1\$CHR[T1\$CHR == "Y"] <- 24
head(T1)
unique(T1\$CHR)
T1\$CHR <- as.numeric(T1\$CHR)


if(nrow(T1)>0) {
${format}("${prefix}.manhattan_mqc.${format}",width = ${((plot_size as double) *3.0 ) as int}, height = ${plot_size}, unit = "px")
manhattan(T1,main="${title}",sub="${subtitle}" ${args1});
dev.off()

${format}("${prefix}.qqplot_mqc.${format}",width = ${plot_size}, height = ${plot_size}, unit = "px")
qq(T1\$P,main="${title}",sub="${subtitle}"  ${args2});
dev.off()
}
__EOF__

R --vanilla < TMP/jeter.R


cat << EOF > versions.yml
${task.process}:
    R: \$(R --version | awk '(NR==1) {print \$3;}')
EOF
"""

stub:
    def prefix = task.ext.prefix?:"${meta.id}"
    def format = task.ext.format?:"png"
"""
touch versions.yml ${prefix}.qqplot_mqc.${format} ${prefix}.manhattan_mqc.${format}
"""
}
