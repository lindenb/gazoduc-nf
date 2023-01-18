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

include {moduleLoad;parseBoolean} from '../utils/functions.nf'

process UCSC_LIFTOVER_01 {
tag "${row.bed} ${row.chain}"
input:
	val(meta)
	val(row)
output:
	tuple row,
		path("remapped.bed${parseBoolean(row.compressed?:"true")?".gz":""}"),
		path("unmapped.bed${parseBoolean(row.compressed?:"true")?".gz":""}"),emit:output
	path("version.xml"), emit:version
script:
	def bed = row.bed
	def chain = row.chain
	"""
	hostname 1>&2
	${moduleLoad("ucsc htslib")}
	set -o pipefail

	mkdir TMP
	liftOver "${bed}" ${chain}" remapped.bed unmapped.bed

	if [ ! -z "${parseBoolean(row.compressed?:"true")?"Y":""}" ] ; then

		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n remapped.bed |\
			bgzip > remapped.bed.gz

		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n unapped.bed |\
			bgzip > unapped.bed.gz

		rm remapped.bed unapped.bed
		tabix -p bed remapped.bed.gz
		tabix -p bed unapped.bed.gz
	fi

	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">liftover bed</entry>
		<entry key="bed">${bed}</entry>
		<entry key="chain">${chain}</entry>
	</properties>
	EOF
	"""
	stub:
	"""
	touch unmapped.bed remapped.bed
	touch unmapped.bed.gz remapped.bed.gz
       	echo "<properties/>" > version.xml
	"""
	}
