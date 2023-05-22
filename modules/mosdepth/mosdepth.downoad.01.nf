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

def gazoduc = gazoduc.Gazoduc.getInstance()

gazoduc.make("mosdepth_version","0.3.3").
        description("mosdepth version").
        put()




process MOSDEPTH_DOWNLOAD_01 {
input:
	val(meta)
output:
	path("mosdepth"),emit:executable
	path("version.xml"),emit:version
script:
	def mosdepth_version = meta.mosdepth_version?:"0.3.3"
	def url = "https://github.com/brentp/mosdepth/releases/download/v${mosdepth_version}/mosdepth"
"""
hostname 1>&2
wget -O mosdepth "${url}"
chmod +x mosdepth

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download mosdepth</entry>
        <entry key="mosdepth.version">${mosdepth_version}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
stub:
"""
touch mosdepth
echo "<properties/>" > version.xml
"""
}
