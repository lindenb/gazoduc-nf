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

include {getVersionCmd} from '../utils/functions.nf'

process GET_LIFTOVER_CHAINS_02 {
tag "${hgFrom} -> ${hgTo}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(hgFrom)
	val(hgTo)
output:
	path("*.chain"),emit:chain
	tuple val(hgFrom),val(hgTo),path("*.chain"),emit:output /* metadata and chain */
	path("version.xml"),emit:version
script:
	def chainName = "${hgFrom}To${hgTo.substring(0,1).toUpperCase()}${hgTo.substring(1)}.over.chain"
	def chain_url = "https://hgdownload.cse.ucsc.edu/goldenpath/${hgFrom}/liftOver/${chainName}.gz"
"""
hostname 1>&2
mkdir -p TMP

wget  -O "TMP/${chainName}.gz" "${chain_url}"
gunzip -f "TMP/${chainName}.gz"

mv -v "TMP/${chainName}" ./ 

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download liftover chains</entry>
	<entry key="ref.src">${hgFrom}</entry>
	<entry key="ref.dest">${hgTo}</entry>
	<entry key="chain.url"><url>${chain_url}</url></entry>
	<entry key="versions">${getVersionCmd("wget")}</entry>
</properties>
EOF
"""
}

