include {moduleLoad} from '../utils/functions.nf'

process VALIDATE_CASE_CONTROL_PED_01 {
tag "${pedigree.name}"
executor "local"
input:
	val(meta)
	path(pedigree)
output:
	path("validated.ped"),emit:pedigree
	path("cases.list"),emit:cases_list
	path("controls.list"),emit:controls_list
	path("version.xml"),emit:version
script:
"""
hostname 2>&1
${moduleLoad("datamash")}
#set -o pipefail

# check no space
grep -F " " "${pedigree}"| head -n 1 | cat > jeter.txt
test ! -s jeter.txt

# check num columns
awk -F '\t' '(NF<6)' "${pedigree}" > jeter.txt
cat -n jeter.txt
test ! -s jeter.txt

# check no extra white spaces
tr "\t" "${pedigree}" | awk '{S=\$1;gsub(/[ ]/,"",S); if(S!=\$1) print;}' > jeter.txt
cat -n jeter.txt
test ! -s jeter.txt

# 6th column must be case or control
awk -F '\t' '!(\$6=="case" || \$6=="control" || \$6=="affected" || \$6=="unaffected")' "${pedigree}" > jeter.txt
cat -n jeter.txt
test ! -s jeter.txt

# 5th column must be male or female
awk -F '\t' '!(\$5=="male" || \$5=="female" || \$5=="0" || \$5=="unknown" || \$5=="undefined")' "${pedigree}" > jeter.txt
cat -n jeter.txt
test ! -s jeter.txt


# make list of cases
awk -F '\t' '(\$6=="case" || \$6=="affected") {print \$2;}' "${pedigree}" | sort | uniq > cases.list

# make list of controls
awk -F '\t' '(\$6=="control" || \$6=="unaffected") {print \$2;}' "${pedigree}" | sort | uniq > controls.list

# check no common between controls and cases
comm -12 cases.list controls.list > jeter.txt
test ! -s jeter.txt

cp "${pedigree}" validated.ped
rm jeter.txt




##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">validate pedigree</entry>
	<entry key="pedigree">${pedigree}</entry>
        <entry key="pedigree.count">\$(wc -l < validated.ped)</entry>
        <entry key="cases.count">\$(wc -l < cases.list)</entry>
        <entry key="controls.count">\$(wc -l < controls.list)</entry>
	<entry key="table">
EOF

cut -f5,6 validated.ped | sort -T . -t '\t' -k1,1 | datamash crosstab 1,2 |\
	awk -F '\t' 'BEGIN {printf("<table><thead><caption>Sex/Status</caption></thead><tbody>");} {printf("<tr>");for(i=1;i<=NF;i++) {printf("<%s>%s</%s>",(NR==1 || i==1?"th":"td"),\$i,(NR==1 || i==1?"th":"td"));} printf("</tr>\\n");} END{printf("</tbody></table>\\n");}' >> version.xml


cat << EOF >> version.xml
</entry>
</properties>
EOF
"""
}

