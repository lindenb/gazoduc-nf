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



def ucscID(meta) {
	return meta.ucsc_name?:"undefined"
	}




def chainURL(ref1,ref2) {
	def g1 = ucscID(ref1)
	def g2 = ucscID(ref2)
	//title-ize
	def g2t = g2.substring(0,1).toUpperCase()+g2.substring(1)
	return "https://hgdownload.cse.ucsc.edu/goldenpath/${g1}/liftOver/${g1}To${g2t}.over.chain.gz"
	}


/** 

get UCSC liftover chain to go from row.ref1 to row.ref2

*/
process DOWNLOAD_CHAIN {
label "process_single"
tag "${meta1.id} to ${meta2.id}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),val(dict1)
	tuple val(meta2),val(dict2)
output:
	tuple val(meta1),val(meta2),path("*.chain"),emit:chain
	path("versions.yml"),emit:versions
script:
	def url = task.ext.url?:"${chainURL(meta1,meta2)}"
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"
"""
hostname 1>&2
set -o pipefail
mkdir -p TMP

curl -L "${url}" |\\
	gunzip -c |\\
	jvarkit ${jvm} convertliftoverchain \\
	-R1 "${dict1}" \\
	-R2 "${dict2}" > TMP/liftover.chain

test -s TMP/liftover.chain

mv TMP/liftover.chain  "${meta1.id}_to_${meta2.id}.chain"

cat << EOF > versions.yml
"${task.process}":
	url:"${url}"
EOF
"""

stub:
"""
touch versions.yml "${meta1.id}_to_${meta2.id}.chain"
"""
}

