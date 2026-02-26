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

process AVERAGE_PIHAT {
tag "${genome.name}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
    tuple val(meta ), path(genome)
    tuple val(meta2), path(sample2group)
output:
    tuple val(meta),path("*.tsv"),optional:true,emit:tsv
    tuple val(meta),path("*.png"),optional:true,emit:png
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"sample2avg.pihat"
    def sub_title = task.ext.sub?:""
    def max_pihat = task.ext.max_pihat?:0.1
"""
mkdir -p TMP

cat << '__EOF__' > TMP/jeter.awk
(NR>1) {
    P[\$2]+=1.0*(\$10);
    P[\$4]+=1.0*(\$10);
    C[\$2]++;
    C[\$4]++;
    }
END {
    for(S in P) {
        printf("%s\t%f\\n",S,P[S]/C[S]);
        }
    }
__EOF__

# create table sample/avg(pihat)/status
awk -f TMP/jeter.awk '${genome}' |\\
	LC_ALL=C sort -T . -t '\t' -k2,2gr > "TMP/jeter.tsv"


cat << '__EOF__' > TMP/jeter.R
T1<-read.table("TMP/jeter.tsv",sep="\\t",header=FALSE,col.names=c("S","X"),colClasses=c("character","numeric"))

# Read the sample-to-group mapping
sample2group <- read.table("${sample2group}", sep="\t", header=FALSE,col.names=c("S","G"), colClasses=c("character", "character"))

# Merge to get group information for each sample
T1 <- merge(T1, sample2group, by="S", all.x=TRUE)

# Assign a color to each group
groups <- unique(T1\$G)
group_colors <- setNames(rainbow(length(groups)), groups)
T1\$color <- group_colors[T1\$G]

png("TMP/jeter.png", width=800, height=800, units="px")
boxplot(T1\$X ,
    ylim=c(0,max(T1\$X)),
    main="AVG(PIHAT)/SAMPLE",
    sub="${sub_title}",
    xlab="Sample",
    ylab="pihat",
    col=T1\$color
    )
abline(h=${max_pihat},col="blue")


# Add legend
legend("topright", legend=groups, fill=group_colors, title="POP")


dev.off()

__EOF__

R --no-save < TMP/jeter.R

mv TMP/jeter.tsv ${prefix}.tsv
mv TMP/jeter.png ${prefix}.png || true

cat << EOF > versions.yml
${task.process}:
    R: todo
EOF
"""
}