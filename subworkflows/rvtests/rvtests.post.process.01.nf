include {R_INSTALL_PACKAGES_01} from '../../modules/R/install.packages.01.nf'
include {getKeyValue;getModules} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'


/** extract test name from filename */
String extractTestName(f) {
    String s = f.getName();
    if(s.endsWith(".assoc")) s=s.substring(0,s.length()-6);
    int dot  = s.lastIndexOf('.');
    s = s.substring(dot+1);
    return s;
    }


workflow RVTESTS_POST_PROCESS {
	take:
		meta
		reference
		subtitle
		list
	main:
		version_ch = Channel.empty()
		R_ch =  R_INSTALL_PACKAGES_01(meta.plus("R_packages":"\"qqman\""))
		version_ch = version_ch.mix(R_ch.version)
		

		assoc2file_ch = list.splitCsv(header: false,sep:'\t',strip:true).
			map{T->T[0]}.
			map{F->[extractTestName(file(F)),F]}.
			groupTuple()
		

		group_ch = groupByTest(meta, assoc2file_ch)
		version_ch = version_ch.mix(group_ch.version)

		plot_ch = plotIt(meta, R_ch.lib, subtitle, group_ch.assoc)
		version_ch = version_ch.mix(plot_ch.version)

		zip_ch = zipIt(meta, plot_ch.plots.concat(group_ch.assoc).collect())

                version_ch = MERGE_VERSION(meta, "rvtests","${subtitle}", version_ch.collect())
	emit:
		zip = zip_ch.zip
		version = version_ch.version
	}
                        


process groupByTest {
tag "${assoc} N=${L.size()}"
afterScript "rm -f jeter.tsv"
input:
	val(meta)
	tuple val(assoc),val(L)
output:
	tuple val(assoc),path("${prefix}${assoc}.tsv"),emit:assoc
	path("version.xml"),emit:version
script:
	prefix = getKeyValue(meta,"prefix","")
"""
hostname 1>&2

cat << EOF > jeter.list
${L.join("\n")}
EOF

# get best sorting key
KEY=`head -n1 '${L[0]}' | awk -F '\t' 'BEGIN{K=""} (NR==1) {for(j=1;j<=3 && K=="";j++) {for(i=1;i<=NF && K=="";i++) {if((\$i=="Pvalue" && j==1) || (\$i=="PermPvalue" && j==2) || (\$i=="PvalueTwoSide" && j==3)){K=sprintf("-k%d,%dg",i,i);break;}}}} END {print K}' `

echo "Sort key will be \${KEY}" 1>&2

head -n 1 '${L[0]}' > jeter.tsv

xargs -a jeter.list -L 50 cat |\
    grep -v -E  '^(Gene|CHROM|Range)\t' |\
    LC_ALL=C sort -t '\t' -T . \${KEY} | uniq >> jeter.tsv

mv jeter.tsv "${prefix}${assoc}.tsv"

cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="Analysis">${assoc}</entry>
  <entry key="Number of files">${L.size()}</entry>
  <entry key="Output">${prefix}${assoc}.tsv</entry>
</properties>
EOF

rm jeter.list
"""
}


process plotIt {
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
	prefix = getKeyValue(meta,"prefix","")
"""
module load ${getModules("r/3.6.3")}
hostname 1>&2
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

############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="name">${task.process}</entry>
  <entry key="Analysis"><dd>${assoc}</entry>
</properties>
EOF


find \${PWD} -type f -name "*.png" >> paths.txt
"""
}

process zipIt {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("${prefix}rvtest.zip"),emit:zip
script:
	prefix = getKeyValue(meta,"prefix","")
"""
cat ${L.join(" ")} | zip -j -9 -@ "${prefix}rvtest.zip"
"""
}
