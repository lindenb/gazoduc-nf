

process RVTESTS_GROUP_BY_TEST_01 {
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

xargs -a jeter.list -L 1 cat |\
    grep -v -E  '^(Gene|CHROM|Range)\t' |\
    LC_ALL=C sort -t '\t' -T . \${KEY} | uniq >> jeter.tsv

mv jeter.tsv "${prefix}${assoc}.tsv"

cat << EOF > version.xml
<properties id="${task.process}">
  <entry key="id">${task.process}</entry>
  <entry key="analysis">${assoc}</entry>
  <entry key="number of files">${L.size()}</entry>
  <entry key="output">${prefix}${assoc}.tsv</entry>
</properties>
EOF
"""
stub:
"""
touch "${prefix}${assoc}.tsv"
echo "<properties/>" > version.xlm
"""
}
