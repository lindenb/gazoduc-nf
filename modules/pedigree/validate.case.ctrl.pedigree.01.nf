/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {moduleLoad} from '../utils/functions.nf'

process VALIDATE_CASE_CONTROL_PED_01 {
tag "${pedigree.name}"
executor "local"
input:
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
awk -F '\t' '{S=\$1;gsub(/[ ]/,"",S); if(S!=\$1) print;}' "${pedigree}" > jeter.txt
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

