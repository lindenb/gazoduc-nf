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

include {moduleLoad;parseBoolean;assertNotEmpty;getKeyValue} from '../utils/functions.nf'

process BED_CLUSTER_01 {
	tag "${bed.name} method:${meta.bed_cluster_method}"
	afterScript "rm -rf TMP"
	memory "2g"
	input:
		val(meta)
		val(genomeId)
		path(bed)
	output:
		path("clusters.list"),emit:output
		path("version.xml"),emit:version
	script:
		if(!meta.containsKey("bed_cluster_method")) throw new IllegalArgumentException("bed_cluster_method undefined");
		def method = meta.bed_cluster_method
		assertNotEmpty(method,"method for bedcluter must be defined")
		def by_chrom = parseBoolean(meta.by_chromosome?:false)
		def chrom_arg = (by_chrom?"--chromosome":"")
		def names_arg = (parseBoolean(meta.with_names?:false)?"--names":"")
		def reference = params.genomes[genomeId].fasta
	"""
	hostname 1>&2
	set -o pipefail
	${moduleLoad("jvarkit")}

	test ! -z "${method}"

	mkdir -p TMP
	java  -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar bedcluster \
		-R "${reference}" \
		${chrom_arg} \
		${method} \
		${names_arg} \
		-o TMP "${bed}"

	mv TMP BEDS

	find \${PWD}/BEDS -type f -name "*.bed" > clusters.list
	test -s clusters.list

	#####################################################################################

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Cluster BED file</entry>
		<entry key="bed">${bed}</entry>
		<entry key="reference">${reference}</entry>
		<entry key="method">${method}</entry>
		<entry key="jvarkit.version">\$(java -jar \${JVARKIT_DIST}/jvarkit.jar --version)</entry>
	</properties>
	EOF
	"""
	stub:
	"""
	touch "clusters.list"
	echo "<properties/>" > version.xml
	"""
	}
