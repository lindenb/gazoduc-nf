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

include {getVersionCmd;isHg19;isHg38;isHs37d5;moduleLoad} from '../utils/functions.nf'

def ucscID(ref) {
	if(isHg19(ref)) return "hg19";
	if(isHg38(ref)) return "hg38";
	throw new IllegalArgumentException("Don't know genome ID at ucsc for "+ref);
	}

def chainURL(ref1,ref2) {
	def g1 = ucscID(ref1)
	def g2 = ucscID(ref2)
	//title-ize
	def g2t = g2.substring(0,1).toUpperCase()+g2.substring(1)
	return "https://hgdownload.cse.ucsc.edu/goldenpath/${g1}/liftOver/${g1}To${g2t}.over.chain.gz"
	}

def renameContig(ref) {
	if(isHs37d5(ref)) return ["https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh37_gencode2UCSC.txt", true ];
	return [];
	}

/** 

get UCSC liftover chain to go from row.ref1 to row.ref2

*/
process GET_LIFTOVER_CHAINS {
tag "${row.reference1} ${row.reference2}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(row)
output:
	path("${row.name?:"liftover"}.chain"),emit:chain /* chain only */
	tuple val(row),path("${row.name?:"liftover"}.chain"),emit:output /* metadata and chain */
	path("version.xml"),emit:version
script:
	def ref1 =  row.reference1
	def ref2 =  row.reference2
	def chain_url = chainURL(ref1,ref2)
	def mapping1 = renameContig(ref1)
	def mapping2 = renameContig(ref2)
"""
hostname 1>&2
${moduleLoad("jvarkit")}
mkdir -p TMP
set -o pipefail

if ${mapping1.size()==2} ; then

	wget -O - "${mapping1[0]}" |\
		awk -F '\t' '{if(\$1=="" || \$2=="") next; printf("%s\t%s\\n",${mapping1[1]?"\$2,\$1":"\$1,\$2"});}' |\
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${ref1}" --column 2 --convert SKIP 	> TMP/src.mapping

	test -s TMP/src.mapping
fi

if ${mapping2.size()==2} ; then

	wget -O - "${mapping2[0]}" |\
		awk -F '\t' '{if(\$1=="" || \$2=="") next; printf("%s\t%s\\n",${mapping2[1]?"\$2,\$1":"\$1,\$2"});}' |\
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${ref2}" --column 2 --convert SKIP 	> TMP/dest.mapping

	test -s TMP/dest.mapping

fi

wget -q -O - "${chain_url}" |\
	gunzip -c |\
	java -jar \${JVARKIT_DIST}/convertliftoverchain.jar \
	-R1 "${mapping1.size()==2?"TMP/src.mapping":ref1}" \
	-R2 "${mapping2.size()==2?"TMP/dest.mapping":ref2}" > TMP/liftover.chain

test -s TMP/liftover.chain

mv TMP/liftover.chain  "${row.name?:"liftover"}.chain"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">Download liftover chains</entry>
	<entry key="ref.src">${ref1}</entry>
	<entry key="ref.dest">${ref2}</entry>
	<entry key="chain.url"><url>${chain_url}</url></entry>
	<entry key="rename1"><url>${mapping1.size()==2?mapping1[0]:""}</url></entry>
	<entry key="rename2"><url>${mapping2.size()==2?mapping2[0]:""}</url></entry>
	<entry key="versions">${getVersionCmd("jvarkit/convertliftoverchain jvarkit/bedrenamechr awk wget")}</entry>
</properties>
EOF
"""
}

