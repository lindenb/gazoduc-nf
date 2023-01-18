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
include {parseBoolean;isBlank} from './functions.nf'

process COLLECT_TO_FILE_01 {
executor "local"
afterScript "rm jeter.list"
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("concat.list"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail

if [ ! -z "${isBlank(meta.header)?"":"Y"}" ] ; then
	echo '${meta.header?:""}' > concat.list
fi

cat << EOF > jeter.list
${L.join("\n")}
EOF

if [ ! -z "${parseBoolean(meta.sort)?"Y":""}" ] ; then

	awk -F '/' '{printf("%s\t%s\\n",\$NF,\$0);}' jeter.list |\
		sort -t '\t' -T. -k1,1 -k2,2 | cut -f 2- | uniq >> concat.list

else

	cat jeter.list >> concat.list

fi


## if it's too fast, prevent clock problems

sleep 10
touch -c concat.list

cat << EOF > version.xml
<properties id='${task.process}'>
  <entry key="name">${task.process}</entry>
  <entry key="description">collecte strings and write it into a file</entry>
  <entry key="count">${L.size()}</entry>
  <entry key="header">${meta.header?:"(no-header)"}</entry>
</properties>
"""
stub:
"""
touch concat.list
echo "<properties/>" > version.xml
"""
}
