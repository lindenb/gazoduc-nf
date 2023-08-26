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

process WGET_01 {
tag "${row.url}"
afterScript "rm -rf TMP"
input:
        val(meta)
        val(row)
output:
        tuple val(meta),path("${row.output}"),emit:output
	path("version.xml"),emit:version
script:
	if(!row.containsKey("url")) throw new IllegalArgumentException("row.url is missing");
	if(!row.containsKey("output")) throw new IllegalArgumentException("row.output is missing");
	def url = row.url;
	def output = row.output;
	if(url.isEmtpy())  throw new IllegalArgumentException("row.url is empty");
	if(output.isEmtpy())  throw new IllegalArgumentException("row.output is empty");

"""
hostname 1>&2

mkdir -p TMP
wget -O 'TMP/${output}' '${url}'
mv -v 'TMP/${output}' ./

#############################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="Name">${task.process}</entry>
	<entry key="description">wget</entry>
        <entry key="url">${url}</entry>
</properties>
EOF

"""
}

