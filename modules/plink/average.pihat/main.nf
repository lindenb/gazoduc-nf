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
tag "${meta.id}"
afterScript "rm -rf TMP"
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
    tuple val(meta ), path(genome) // plink genome-file
    tuple val(meta2), path(sample2group)
output:
    tuple val(meta),path("*.avg.pihat.tsv"),emit:tsv
    tuple val(meta),path("*.exclude_mqc.tsv"),emit:exclude
    tuple val(meta),path("*.png"),optional:true,emit:png
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}"
    def sub_title = task.ext.sub?:""
    def max_pihat = task.ext.max_pihat?:0.1
    def plot_size= task.ext.plot_size?:1000
    def format = task.ext.format?:"png"
"""
mkdir -p TMP

# check genome file is ok
awk '(NR==1) {print \$10;}' "${genome}"  | grep -w -F PI_HAT 1>&2

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
    sort -T TMP -t '\t' -k1,1 > "TMP/jeter.tsv"

if ${sample2group?true:false}
then

    sort -T TMP -t '\t' -k1,1 --unique "${sample2group}" > TMP/jeter.a.tsv
    join -t '\t' -1 1 -2 1 -e other -a 1 -o '1.1,1.2,2.2'  TMP/jeter.tsv TMP/jeter.a.tsv > TMP/jeter2.tsv
    mv TMP/jeter2.tsv  TMP/jeter.tsv

else

    awk -F '\t' '{printf("%s\tall\\n",\$0);}' TMP/jeter.tsv > TMP/jeter2.tsv
    mv TMP/jeter2.tsv TMP/jeter.tsv
fi

sort -T TMP -t '\t' -k2,2gr  "TMP/jeter.tsv"  > TMP/jeter2.tsv
mv TMP/jeter2.tsv TMP/jeter.tsv

head TMP/jeter.tsv 1>&2

cat << '__EOF__' > TMP/jeter.R
T1<-read.table("TMP/jeter.tsv",sep="\\t",header=FALSE,stringsAsFactors=FALSE, col.names=c("sample_name","p_value","group_name"),colClasses=c("character","numeric","character"))
head(T1)

# create a vector or unique groups
groups <- unique(T1\$group_name)

cols = rainbow(length(groups))
${format}("TMP/jeter.${format}",width = ${plot_size}, height = ${plot_size}, unit = "px")


# Draw boxplot by group
bp <- boxplot(p_value ~ group_name,
    data = T1,
    col = cols,
    border = "gray30",
    main = "Average PIHAT",
    ylab = "average(p_value)",
    xlab = "Group",
    las = 2)

if (length(groups) > 1) {
  legend("topright",
         legend = groups,
         fill = cols,
         border = "gray30",
         bty = "n",
         cex = 0.9)
}


dev.off()

__EOF__

R --no-save < TMP/jeter.R

awk -F '\t' '(\$2*1.0 > ${max_pihat})' TMP/jeter.tsv > ${prefix}.exclude_mqc.tsv
mv TMP/jeter.tsv ${prefix}.avg.pihat.tsv
mv TMP/jeter.${format} ${prefix}.avg.pihat_mqc.${format} || true

cat << EOF > versions.yml
${task.process}:
    R: \$(R --version | awk '(NR==1) {print \$3;}')
EOF
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}"
"""
touch versions.yml ${prefix}.tsv ${prefix}.png
"""
}
