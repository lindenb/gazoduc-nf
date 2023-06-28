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

log.info("parsing picard.scatter.nf ....")

include {moduleLoad} from '../utils/functions.nf'


process SCATTER_INTERVALS_BY_NS {
tag "${fasta.name}"
memory "3g"
input:
	val(meta)
	path(fasta)
output:
	path("${fasta.getSimpleName()}.${meta.OUTPUT_TYPE}.${meta.MAX_TO_MERGE}.interval_list"), emit:output
	path("version.xml"), emit:version

script:
	if(!meta.containsKey("OUTPUT_TYPE")) throw new IllegalArgumentException("meta.OUTPUT_TYPE is missing");
	if(!meta.containsKey("MAX_TO_MERGE")) throw new IllegalArgumentException("meta.MAX_TO_MERGE is missing");


	type = meta.OUTPUT_TYPE
	maxToMerge = meta.MAX_TO_MERGE

	"""
	hostname 1>&2
	${moduleLoad("picard")}

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=. -jar \${PICARD_JAR} ScatterIntervalsByNs \
	    -R "${fasta.toRealPath()}" \
	    --MAX_TO_MERGE "${maxToMerge}" \
	    -O "${fasta.getSimpleName()}.${type}.${maxToMerge}.interval_list" \
	    -OUTPUT_TYPE "${type}"

	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">call ScatterIntervalsByNs to get a list of intervals</entry>
		<entry key="reference">${fasta}</entry>
		<entry key="type">${type}</entry>
		<entry key="maxToMerge">${maxToMerge}</entry>
		<entry key="picard">\$(java -jar \${PICARD_JAR} ScatterIntervalsByNs --version 2>&1)</entry>
	</properties>
	EOF
	"""
	}

