/*

Copyright (c) 2024 Pierre Lindenbaum

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

process MERGE_VERSION {
tag "N=${L.size()}"
executor "local"
input:
	val(name)
	val(L)
output:
	path("${params.prefix?:""}version.xml"),emit:version
script:
	prefix = params.prefix?:""
"""
cat << EOF > jeter.xml
<properties id="${name}">
	<entry key="name">${name}</entry>
	<entry key="date">\$(date)</entry>
	<entry key="steps">
EOF

for X in ${L.join(" ")}
do
	xmllint --nocdata --format "\${X}" | tail -n+2 >> jeter.xml
done

cat << EOF >> jeter.xml
	</entry>
</properties>
EOF
	xmllint --format jeter.xml > "${prefix}version.xml"
rm jeter.xml
"""
stub:
"""
echo "<properties/>" > "${prefix}version.xml"
"""
}
