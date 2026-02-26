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
include {verify            } from '../../../modules/utils/functions.nf'
include {isBlank           } from '../../../modules/utils/functions.nf'

process PLOT_MDS {
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta1),path(sample2pop) //optional, TSV
    tuple val(meta),path(mds)
output:
    tuple val(meta),path("*.png"),optional:true,emit:png
    tuple val(meta),path("*.pdf"),optional:true,emit:pdf
    path("versions.yml"),emit:versions
script:
    def format= task.ext.format?:"png"
    verify((format=="png" || format=="pdf"),"${task.process} format should be png or pdf but got ${format}")
    def Cx = task.Cx?:"${meta.Cx}"
    verify(!isBlank(Cx),"${task.process} Cx shouldn't be blank")
    def Cy = task.Cy?:"${meta.Cy}"
    verify(!isBlank(Cy),"${task.process} Cy shouldn't be blank")
    def prefix= task.ext.prefix?:"${meta.id}.mds"
    def R_main = task.ext.R_main?:"PCA"
    def R_sub = task.ext.R_sub?:"${Cx} / ${Cy}"
    def plot_size = task.ext.plot_size?:1000
"""
mkdir -p TMP
cat << '__EOF__' > TMP/jeter.R

mds_file <- "${mds}"
pop_file <- "${sample2pop?:"NO_SAMPLE_TO_POP"}"

mds <- read.table(mds_file,header=TRUE, stringsAsFactors = FALSE)
mds\$color <-rep("black", nrow(mds))


## If population file exists, read + merge population info
if (file.exists(pop_file)) {
    # Read the sample-to-group mapping
    sample2group <- read.table(pop_file, sep="\t", header=FALSE,col.names=c("IID","G"), colClasses=c("character", "character"), stringsAsFactors = FALSE)

    # Merge to get group information for each sample
    mds <- merge(mds, sample2group, by="IID", all.x=TRUE, sort = FALSE)

    # Assign a color to each group
    groups <- unique(sample2group\$G)
    group_colors <- setNames(rainbow(length(groups)), groups)
    mds\$color <- group_colors[mds\$G]
}

head(mds)

${format}("TMP/jeter.${format}",width = ${plot_size}, height = ${plot_size}, unit = "px")
plot(
    mds\$${Cx},
    mds\$${Cy},
    main ="${R_main}",
    pch = 19,
    xlab = "${Cx}",
    ylab = "${Cy}",
    sub="${R_sub}",
    col=mds\$color
    )


## Legend only if we actually have populations/colors
if (file.exists(pop_file)) {
    # Add legend
    legend("topright", legend=groups, fill=group_colors, title="POP")
}

dev.off()
__EOF__






R --no-save < TMP/jeter.R

mv TMP/jeter.${format} "${prefix}.${Cx}_${Cy}_mqc.${format}"

cat << EOF > versions.yml
${task.process}:
    R: \$(R --version | awk '(NR==1) {print \$3;}')
EOF

"""
stub:
    def format= task.ext.format?:"png"
    verify((format=="png" || format=="pdf"),"${task.process} format should be png or pdf but got ${format}")
    def Cx = task.Cx?:"${meta.Cx}"
    verify(!isBlank(Cx),"${task.process} Cx shouldn't be blank")
    def Cy = task.Cy?:"${meta.Cy}"
    verify(!isBlank(Cy),"${task.process} Cy shouldn't be blank")
    def prefix= task.ext.prefix?:"${meta.id}.mds"
"""
touch versions.yml "${prefix}.${Cx}_${Cy}_mqc.pdf"  "${prefix}.${Cx}_${Cy}_mqc.png"
"""
}