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

process PEDIGREE_FOR_RVTESTS {
tag "${pedigree.name}"
executor "local"
input:
        val(meta)
        path(pedigree)
output:
        path("rvtests.pedigree.ped"),emit:pedigree
        path("version.xml"),emit:version
script:
	def p= params.rvtests.phenotype_name
"""
hostname 1>&2
set -o pipefail

# assert no duplicate between cases and controls
cut -f 2 "${pedigree}" | sort | uniq -d  > dups.txt
test ! -s dups.txt
rm dups.txt

awk -F '\t' 'BEGIN {print "fid\tiid\tfatid\tmatid\tsex\t${p}"} (\$6=="case" || \$6=="affected") {printf("%s\t%s\t0\t0\t0\t2\\n",\$1,\$2);next;} (\$6=="control" || \$6=="unaffected") {printf("%s\t%s\t0\t0\t0\t1\\n",\$1,\$2);next;}' "${pedigree}" > jeter.ped


############################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">create pedigree for rvtests</entry>
        <entry key="pedigree">${pedigree.toRealPath()}</entry>
        <entry key="cases.count">\$(awk -F '\t' '(\$6==2)' jeter.ped | wc -l)' )</entry>
        <entry key="controls.count">\$(awk -F '\t' '(\$6==1)' jeter.ped | wc -l)' )</entry>
</properties>
###############################################
EOF

mv -v jeter.ped rvtests.pedigree.ped
"""
}
