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

include {moduleLoad;isUrl;isBlank;getModules;parseBoolean;getGencodeGff3Url} from '../utils/functions.nf'

process DOWNLOAD_GFF3_01 {
tag "${file(reference).name}"
afterScript "rm -rf TMP"
memory "1g"
input:
	val(meta)
	val(reference)
output:
	path("${file(reference).getSimpleName()}.gff3${parseBoolean(meta.with_tabix?:"true")?".gz":""}"),emit:gff3
	path("${file(reference).getSimpleName()}.gff3.gz.tbi"),emit:tbi,optional:true
	path("version.xml"),emit:version
script:
	def url0 = meta.gff3url?:""
	def url = (isBlank(url0)?getGencodeGff3Url(reference):url0)
	def with_tabix = parseBoolean(meta.with_tabix?:"true")

	"""
	hostname 1>&2
	${moduleLoad("htslib jvarkit")}
	set -o pipefail
	mkdir TMP

	if [ ! -z "${url.equals("undefined") || isBlank(url)?"Y":""}" ] ; then
		touch TMP/jeter0.gff3
	elif [ ! -z "${isUrl(url)?"Y":""}" ] ; then
		wget -O TMP/jeter0.gff3 "${url}"
	else
		cp -v "${url}" TMP/jeter0.gff3
	fi

	if [[ `file TMP/jeter0.gff3 | grep gzip` ]] ; then
		mv TMP/jeter0.gff3 TMP/jeter0.gff3.gz
	fi
	

	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar --column 1 --convert SKIP -R "${reference}" TMP/jeter0.* > TMP/jeter.gff3
	rm TMP/jeter0.*
	mv TMP/jeter.gff3 TMP/jeter0.gff3


	
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -k1,1 -k4,4n TMP/jeter0.gff3 > TMP/jeter.gff3
	rm TMP/jeter0.gff3

	if [ ! -z "${with_tabix?"Y":""}" ] ; then
		bgzip TMP/jeter.gff3
		tabix -p gff TMP/jeter.gff3.gz
		mv TMP/jeter.gff3.gz "${file(reference).getSimpleName()}.gff3.gz"
		mv TMP/jeter.gff3.gz.tbi "${file(reference).getSimpleName()}.gff3.gz.tbi"
	else
		mv TMP/jeter.gff3 "${file(reference).getSimpleName()}.gff3"		
	fi
	
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download GFF for Reference</entry>
		<entry key="reference">${reference}</entry>
		<entry key="with_tabix">${with_tabix}</entry>
		<entry key="url"><a>${url}</a></entry>
		<entry key="bedrenamechr">\$(java  -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
		<entry key="bgzip">\$(bgzip --version | head -n 1)</entry>
		<entry key="tabix">\$(tabix --version | head -n 1)</entry>
	</properties>
	EOF
	"""
	
	stub:
	"""
	touch "${file(reference).getSimpleName()}.gff3${parseBoolean(meta.with_tabix?:"true")?".gz":""}"
        touch "${file(reference).getSimpleName()}.gff3.gz.tbi"
	echo "<properties/>" > version.xml
	"""
	}
