process RVTESTS_PLOT_ASSOC_01 {
tag "${assoc} (${tsv.name})"
afterScript "rm -f jeter.tsv jeter.R"
input:
	val(meta)
	path(rlib)
	val(subtitle)
	tuple val(assoc),val(tsv)
output:
	path("paths.txt"),emit:plots
	path("version.xml"),emit:version
script:
	def r_version= getKeyValue(meta,"R_module","r/3.6.3")
	prefix = getKeyValue(meta,"prefix","")
"""
hostname 1>&2
${moduleLoad(r_version)}
realpath "${tsv}" > paths.txt

awk -F '\t' 'BEGIN{K=-1; R=-1;FS="\t";} (NR==1) {for(i=1;i<=NF && R<0;i++) {if(\$i=="RANGE") {R=i;}}for(j=1;j<=3 && K<0;j++) {for(i=1;i<=NF && K<0;i++) {if((\$i=="Pvalue" && j==1) || (\$i=="PermPvalue" && j==2) || (\$i=="PvalueTwoSide" && j==3)){K=i;break;}}} printf("SNP\tCHR\tBP\tP\\n");next;} {if(K<0 || R<0 || \$K=="NA") next;split(\$R,a,/[:-]/);C=a[1];gsub("^chr","",C);if(C=="X") {C="23";} else if(C=="Y") {C="24";} if(int(C)<0) next; printf("%s\t%s\t%s\t%s\\n",\$1,C,a[2],\$K);}' '${tsv}' > jeter.tsv

wc -l jeter.tsv
head jeter.tsv

cat << '__EOF__' > jeter.R
library("qqman",lib.loc="${rlib.toRealPath()}")
T1 <- read.table("jeter.tsv",header=TRUE,sep="\t",stringsAsFactors=FALSE)

if(nrow(T1)>0) {
png("${prefix}${assoc}.manhattan.png")
manhattan(T1,main="${prefix}${assoc}",sub="${subtitle}");
dev.off()

png("${prefix}${assoc}.qqplot.png")
qq(T1\$P,main="${prefix}${assoc}",sub="${subtitle}");
dev.off()
}
__EOF__

R --vanilla < jeter.R || true

cat << EOF > version.xml
<dl id="${task.process}">
  <dt>name</dt><dd>${task.process}</dd>
  <dt>Analysis</dt><dd>${assoc}</dd>
</dl>
EOF


find \${PWD} -type f -name "*.png" >> paths.txt

#########################################################################
cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="id">${task.process}</entry>
  <entry key="description">plot output of rvtests</entry>
  <entry key="analysis">${assoc}</entry>
  <entry key="R.version">\$(R --version | head -n1)</entry>
</properties>
EOF
"""
stub:
"""
touch paths.txt
echo "<properties/>" > version.xml
"""
}
