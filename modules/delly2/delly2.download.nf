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

include {getVersionCmd} from '../utils/functions.nf'


process DOWNLOAD_DELLY2 {
	errorStrategy "retry"
	maxRetries 3
	output:
		path("delly"),emit:output
		path("version.xml"),emit:version
	script:
		def version = params.delly2.version
		def url = "https://github.com/dellytools/delly/releases/download/v${version}/delly_v${version}_linux_x86_64bit"
	"""
	hostname 1>&2
	wget -O delly "${url}"
	touch -c delly
	chmod a+x delly
	sleep 10

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download delly</entry>
		<entry key="version">${version}</entry>
		<entry key="url"><a>${url}</a></entry>
		<entry key="version">${getVersionCmd("wget")}</entry>
	</properties>
	EOF
	"""
	}
