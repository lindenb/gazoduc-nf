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
process RVTESTS_PLOT_ASSOC {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
	tuple val(meta),val(tsv)
output:
	tuple val(meta),path("*.manhattan_mqc.png"),optional:true,emit:manhattan
	tuple val(meta),path("*.qqplot_mqc.png"),optional:true,emit:qqplot
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
	def title = task.ext.title?:"${meta.id}"
	def subtitle = task.ext.subtitle?:""
"""
hostname 1>&2

awk -F '\t' 'BEGIN{K=-1; R=-1;FS="\t";} (NR==1) {for(i=1;i<=NF && R<0;i++) {if(\$i=="RANGE") {R=i;}}for(j=1;j<=3 && K<0;j++) {for(i=1;i<=NF && K<0;i++) {if((\$i=="Pvalue" && j==1) || (\$i=="PermPvalue" && j==2) || (\$i=="PvalueTwoSide" && j==3)){K=i;break;}}} printf("SNP\tCHR\tBP\tP\\n");next;} {if(K<0 || R<0 || \$K=="NA") next;split(\$R,a,/[:-]/);C=a[1];gsub("^chr","",C);if(C=="X") {C="23";} else if(C=="Y") {C="24";} if(int(C)<0) next; printf("%s\t%s\t%s\t%s\\n",\$1,C,a[2],\$K);}' '${tsv}' > jeter.tsv

wc -l jeter.tsv
head jeter.tsv

cat << '__EOF__' > jeter.R
library("qqman")
T1 <- read.table("jeter.tsv",header=TRUE,sep="\t",stringsAsFactors=FALSE)

if(nrow(T1)>0) {
png("${prefix}.manhattan_mqc.png")
manhattan(T1,main="${title}",sub="${subtitle}");
dev.off()

png("${prefix}.qqplot_mqc.png")
qq(T1\$P,main="${title}",sub="${subtitle}");
dev.off()
}
__EOF__

R --vanilla < jeter.R || true

cat << EOF > versions.yml
${task.process}:
	R: \$(R --version | head -n1)
EOF
"""
stub:
"""
touch versions.yml
"""
}
