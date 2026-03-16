/*
Copyright (c) 2026 Pierre Lindenbaum

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
include { isBlank } from "../../../modules/utils/functions.nf"
include { verify  } from "../../../modules/utils/functions.nf"


process DOWNLOAD_GENCODE {
label "process_single"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(dict)
output:
	tuple val(meta),path("*.txt.gz"),path("*.txt.gz.tbi"),path("*.sql"),emit:gencode
	path("versions.yml"),emit:versions
script:
	def prefix=task.ext.prefix?:"${meta.id}.gencode"
	def version= task.ext.version?:"V49"
	def jvm = task.ext.jvm?:"-Djava.io.tmpdir=TMP "
	def url= task.ext.url?:""
	if(isBlank(url)) {
		if(meta.ucsc_name=="hg19") {
			url ="https://hgdownload.soe.ucsc.edu/goldenpath/hg19/database/wgEncodeGencodeComp${version}lift37.txt.gz"
			}
		else if(meta.ucsc_name=="hg38") {
			url ="https://hgdownload.soe.ucsc.edu/goldenpath/hg38/database/wgEncodeGencodeComp${version}.txt.gz"
			}
		}
	verify(!isBlank(url),"${task.process} url is empty")
	verify(url.endsWith(".txt.gz") , "${task.process} url ${url} should end with .txt.gz")
"""
hostname 1>&2
mkdir -p TMP

curl -L -o TMP/gencode.txt.gz  "${url}"
curl -L -o TMP/gencode.sql  "${url.replaceAll("\\.txt\\.gz\$",".sql")}"

gunzip -c TMP/gencode.txt.gz |\\
		jvarkit ${jvm} bedrenamechr -f "${dict}" --column 3 --convert SKIP |\\
		sort -S ${task.memory.kilo} -T TMP -t '\t' -k3,3 -k5,5n |\\
		bgzip > "${prefix}.gencode.txt.gz"

tabix -f -0 -b 5 -e 6 -s 3 "${prefix}.gencode.txt.gz"

mv TMP/gencode.sql  "${prefix}.gencode.sql"


cat << END_VERSIONS > versions.yml
"${task.process}":
	url: "${url}"
	version : "${version}"
END_VERSIONS
"""
stub:
def prefix=task.ext.prefix?:"${meta.id}.gencode"
"""
touch versions.yml ${prefix}.txt.gz ${prefix}.txt.gz.tbi 
"""
}

