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
nextflow.enable.dsl=2

include {getKeyValue;getModules} from '../utils/functions.nf'


process INTERVAL_LIST_TO_BED {
tag "${file(interval_list).name}"
label "process_low"
input:
	val(meta)
	val(interval_list)
output:
	path("${file(interval_list).getSimpleName()}.bed"), emit:bed
	path("version.xml"), emit:version
script:
	def SORT = getKeyValue(meta,"SORT","true")
	def SCORE = getKeyValue(meta,"SCORE","500")
"""
	hostname 1>&2
	module load ${getModules("picard")}

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${PICARD_JAR} IntervalListToBed \
		--INPUT "${interval_list}" \
		--OUTPUT "${file(interval_list).getSimpleName()}.bed"  \
		--SCORE ${SCORE} \
		--SORT ${SORT}
		
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">call IntervalListToBed to convert intervals to bed</entry>
	<entry key="interval_list">${interval_list}</entry>
	<entry key="SCORE">${SCORE}</entry>
	<entry key="SORT">${SORT}</entry>
	<entry key="picard">\$(java -jar \${PICARD_JAR} IntervalListToBed --version 2>&1)</entry>
</properties>
EOF
"""
stub:
"""
touch "${file(interval_list).getSimpleName()}.bed"
echo "<properties/>" > version.xml
"""
}
