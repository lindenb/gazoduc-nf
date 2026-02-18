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

process PLOT_MDS {
tag "${mds.name} ${C1} ${C2} ${format}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta),path(mds),val(C1),val(C2),val(format)
output:
    tuple val(meta),path("*.{pdf,png}"),emit:plot
    path("versions.yml"),emit:versions
script:
    def prefix= task.ext.prefix?:meta.id+".mds"
    def R_main = task.ext.R_main?:"PCA"
    def R_sub = "${C1} / ${C2}"
"""
mkdir -p TMP
cat << '__EOF__' > TMP/jeter.R
strat <- read.table(file="${mds}",header=TRUE)

${format}("TMP/jeter.${format}")
plot(strat\$${C1},strat\$${C2},main ="${R_main}",sub="${R_sub}")
dev.off()
__EOF__

R --no-save < TMP/jeter.R

mv TMP/jeter.${format} "${prefix}.${C1}.${C2}.${format}"

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
"""
}