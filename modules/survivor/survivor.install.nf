/*

Copyright (c) 2022 Pierre Lindenbaum

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

process SURVIVOR_INSTALL {
input:
	val(meta)
output:
	path("SURVIVOR/Debug/SURVIVOR"),emit:executable
	path("version.xml"),emit:version
script:
	def url="https://github.com/fritzsedlazeck/SURVIVOR.git"
"""
hostname 1>&2
git clone "${url}"
(cd SURVIVOR/Debug && make)


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download and install survivor</entry>
        <entry key="url"><a>${url}</a></entry>
        <entry key="survivor.version">\$(./SURVIVOR/Debug/SURVIVOR  2>&1 /dev/null | grep '^Version' | cut -d ':' -f 2- )</entry>
</properties>
EOF

"""
stub:
"""
mkdir -p SURVIVOR/Debug
touch SURVIVOR/Debug/SURVIVOR
echo "<properties/>" > version.xml
"""
}
