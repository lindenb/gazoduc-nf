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

include {isUrl;moduleLoad;getBoolean;parseBoolean;getKeyValue;getModules;getGencodeGtfUrl} from '../utils/functions.nf'

process DOWNLOAD_GTF_01 {
tag "${file(reference).name}"
afterScript "rm -rf TMP"
memory "1g"
input:
	val(meta)
	val(reference)
output:
	path("${file(reference).getSimpleName()}.gtf${parseBoolean(getKeyValue(meta,"with_tabix",true))?".gz":""}"),emit:gtf
	path("${file(reference).getSimpleName()}.gtf.gz.tbi"),emit:tbi,optional:true
	path("version.xml"),emit:version
script:

	def url0 = getKeyValue(meta,"gtfurl","")
	def url = url0.isEmpty()?getGencodeGtfUrl(reference):url0
	def with_tabix = parseBoolean(getKeyValue(meta,"with_tabix",true))
	"""
	hostname 1>&2
	${moduleLoad("htslib jvarkit")}
	set -o pipefail
	mkdir TMP

	test ! -z "${url}"
	
	if [[ ! -z "${isUrl(url)?"Y":""}" ]] ;  then
		wget -O TMP/jeter0.gtf "${url}"
	else
		test -s "${url}"
		cp -v "${url}" TMP/jeter0.gtf
	fi

	if [[ `file TMP/jeter0.gtf | grep gzip` ]] ; then
		mv TMP/jeter0.gtf TMP/jeter0.gtf.gz
	fi
	

	java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedrenamechr.jar --column 1 --convert SKIP -R "${reference}" TMP/jeter0.* > TMP/jeter.gtf
	rm TMP/jeter0.*
	mv TMP/jeter.gtf TMP/jeter0.gtf


	
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -k1,1 -k4,4n TMP/jeter0.gtf > TMP/jeter.gtf
	rm TMP/jeter0.gtf

	if [ ! -z "${with_tabix?"Y":""}" ] ; then

		bgzip TMP/jeter.gtf
		tabix -p gff TMP/jeter.gtf.gz
		mv TMP/jeter.gtf.gz "${file(reference).getSimpleName()}.gtf.gz"
		mv TMP/jeter.gtf.gz.tbi "${file(reference).getSimpleName()}.gtf.gz.tbi"
	
	else
		mv TMP/jeter.gtf "${file(reference).getSimpleName()}.gtf"
	fi	

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download GTF for Reference</entry>
		<entry key="reference">${reference}</entry>
		<entry key="url"><a>${url}</a></entry>
		<entry key="bedrenamechr">\$(java  -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
		<entry key="with_tabix">${with_tabix}</entry>
		<entry key="bgzip">\$(bgzip --version | head -n 1)</entry>
		<entry key="tabix">\$(tabix --version | head -n 1)</entry>
	</properties>
	EOF
	"""
stub:
	"""
	touch "${file(reference).getSimpleName()}.gtf${parseBoolean(getKeyValue(meta,"with_tabix",true))?".gz":""})"
	touch "${file(reference).getSimpleName()}.gtf.gz.tbi")"
	echo "<properties/>" > version.xml
	"""
	}
