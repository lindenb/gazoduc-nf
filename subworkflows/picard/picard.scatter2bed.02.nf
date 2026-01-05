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
include {moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'


workflow SCATTER_TO_BED {
	take:
		meta
		rows /* contains reference or genomeId */
	main:
		version_ch = Channel.empty()

		scatter_ch = SCATTER_INTERVALS_BY_NS(meta,rows)
		version_ch = version_ch.mix(scatter_ch.version)

		scatter_ch = scatter_ch.output.map{T->T[0].plus("scatter_interval_list":T[1])}	

		bed_ch = INTERVAL_LIST_TO_BED([:],scatter_ch)
		version_ch = version_ch.mix(bed_ch.version)

		rows = bed_ch.output.map{T->T[0].plus("scatter_bed":T[1])}

		version_ch = MERGE_VERSION("scatterBed", version_ch.collect())
	emit:
		output = rows /* scatter_bed and scatter_interval_list */
		version = version_ch
	}


process SCATTER_INTERVALS_BY_NS {
tag "${row.genomeId?:row.reference}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("*.interval_list"), emit:output
	path("version.xml"), emit:version

script:
	if(row.containsKey("genomeId") && row.containsKey("reference") && !params.genomes[genomeId].equals(row.reference)) throw new IllegalArgumentException("mismatch genomeId reference");
	if(!meta.containsKey("OUTPUT_TYPE")) throw new IllegalArgumentException("meta.OUTPUT_TYPE is missing");
	if(!meta.containsKey("MAX_TO_MERGE")) throw new IllegalArgumentException("meta.MAX_TO_MERGE is missing");

	def fasta
	if(row.containsKey("genomeId")) {
		fasta = file("${params.genomes[row.genomeId].fasta}")
		}
	else  if(row.containsKey("reference")) {
		fasta = file("${row.reference}")
		}
	else {
		throw new IllegalArgumentException("missing: genomeId OR reference");
		}

	"""
	hostname 1>&2
	${moduleLoad("gatk/0.0.0")}
	mkdir -p TMP

	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ScatterIntervalsByNs \
	    -R "${fasta.toRealPath()}" \
	    --MAX_TO_MERGE "${meta.MAX_TO_MERGE}" \
	    -O "TMP/jeter.interval_list" \
	    -OUTPUT_TYPE "${meta.OUTPUT_TYPE}"

	mv -v "TMP/jeter.interval_list" "${fasta.getSimpleName()}.${meta.OUTPUT_TYPE}.${meta.MAX_TO_MERGE}.interval_list"


	cat << EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">call ScatterIntervalsByNs to get a list of intervals</entry>
		<entry key="reference">${fasta}</entry>
	</properties>
	EOF
	"""
	}




process INTERVAL_LIST_TO_BED {
tag "${row.scatter_interval_list.name}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(meta)
	val(row)
output:
	tuple val(row),path("${row.scatter_interval_list.getSimpleName()}.bed"), emit:output
	path("version.xml"), emit:version
script:
	def interval_list = row.scatter_interval_list
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0")}
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" IntervalListToBed \
		--INPUT "${interval_list}" \
		--OUTPUT TMP/jeter.bed  \
		${params.gatk.intervalListToBed.args}

LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n  TMP/jeter.bed > "${interval_list.getSimpleName()}.bed"

		
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">call IntervalListToBed to convert intervals to bed</entry>
	<entry key="interval_list">${interval_list}</entry>
</properties>
EOF
"""
}

