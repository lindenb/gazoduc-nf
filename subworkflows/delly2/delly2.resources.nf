/*

Copyright (c) 2024 Pierre Lindenbaum

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

include { moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include { DOWNLOAD_DELLY2 } from '../../modules/delly2/delly2.download.nf'
include { SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'



workflow DELLY2_RESOURCES {
     take:
        meta /* meta */
       	genomeId
     main:
		version_ch = Channel.empty()
		map_ch = GET_MAPPABILITY([:],genomeId)
		version_ch = version_ch.mix(map_ch.version)


		gaps_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"N","MAX_TO_MERGE":"1"], params.genomes[genomeId].fasta)
		version_ch = version_ch.mix(gaps_ch.version)

		exclude_ch = GET_EXCLUDE([:], genomeId)
		version_ch = version_ch.mix(exclude_ch.version)

		delly2_ch = DOWNLOAD_DELLY2()
		version_ch = version_ch.mix(delly2_ch.version)

		xmerge_ch = MERGE_EXCLUDE([:],genomeId,gaps_ch.bed,exclude_ch.bed)
		version_ch = version_ch.mix(xmerge_ch.version)

	emit:
		executable = delly2_ch.output
		exclude = xmerge_ch.bed
		mappability = map_ch.mappability
		version = version_ch
	}


process GET_MAPPABILITY {
tag "${genomeId}"
input:
	val(meta)
	val(genomeId)
output:
	path("blacklist.gz"),emit:mappability
	path("version.xml"),emit:version
	path("blacklist.gz.fai")
	path("blacklist.gz.gzi")
script:
	def genome = params.genomes[genomeId]
	def url = genome.delly2_blacklist_url?:""
	
	if(!url.isEmpty())
	"""
	wget -O blacklist.gz "${url}"
	wget -O blacklist.gz.fai "${url}.fai"
	wget -O blacklist.gz.gzi "${url}.gzi"

	# be sure there is no clock problem....
	touch -c blacklist.gz
	sleep 5
	touch -c blacklist.gz.fai
	sleep 5
	touch -c blacklist.gz.gzi
	sleep 10

	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download mappability files</entry>
		<entry key="url"><a>${url}</a></entry>
	</properties>
	EOF
	"""
 	else
    	"""
 	echo "Undefined build" 1>&2
 	"""
	}



process GET_EXCLUDE {
tag "${genomeId}"
afterScript "rm -f jeter.bed jeter2.bed jeter.interval_list"
input:
	val(meta)
	val(genomeId)
output:
	path("exclude.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def url1 = genome.delly2_exclude_url?:""
"""
hostname 1>&2
${moduleLoad("jvarkit")}
set -o pipefail

if [ ! -z "${url1}" ] ; then
	wget -O - "${url1}" |\
		cut -f 1,2,3|\
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -R "${reference}" --column 1 --convert SKIP > exclude.bed 

else
		touch exclude.bed
fi

	sleep 10

	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download Exclude File</entry>
		<entry key="url"><a>${url1}</a></entry>
		<entry key="version">${getVersionCmd("wget jvarkit/bedrenamechr")}</entry>
	</properties>
	EOF
"""
}

process MERGE_EXCLUDE {
label "process_tiny"
input:
	val(meta)
	val(genomeId)
	path(xclude)
	path(gaps)
output:
	path("exclude.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
"""
hostname 1>&2
${moduleLoad("bedtools")}

awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {printf("%s\t0\t%s\\n",\$1,\$2);}'  "${reference}.fai" > jeter2.bed

cut -f1-3 ${xclude} ${gaps} jeter2.bed | \
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > exclude.bed

rm jeter2.bed


	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">merge ${xclude} , ${gaps} and non standard contigs from reference</entry>
		<entry key="version">${getVersionCmd("bedtools")}</entry>
	</properties>
	EOF

sleep 10

"""
}
