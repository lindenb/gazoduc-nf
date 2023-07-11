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

process SHAPEIT_DOWNLOAD_TOOLS {
afterScript "rm -f jeter.tar.gz"
input:
	val(meta)
output:
	path("ligate_static"),emit:ligate
	path("phase_common_static"),emit:phase_common
	path("phase_rare_static"),emit:phase_rare
	path("version.xml"),emit:version
script:
	def version = params.shapeit.version
	def url = "https://github.com/odelaneau/shapeit5/releases/download/v${version}"
"""
for S in ligate_static phase_common_static phase_rare_static
do
	wget -O "\${S}" "${url}/\${S}"
	chmod +x "\${S}"
done

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download shapeit</entry>
        <entry key="version">${version}</entry>
</properties>
EOF
"""
}
