/*

Copyright (c) 2022 Pierre Lindenbaum

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

include { getKeyValue; getModules; isHg19; isHg38} from '../../modules/utils/functions.nf'
include { SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'


workflow DELLY2_RESOURCES {
     take:
        meta /* meta */
        reference /* indexed fasta reference */
     main:
		version_ch = Channel.empty()
		map_ch = GET_MAPPABILITY([:],reference)
		version_ch = version_ch.mix(map_ch.version)


		gaps_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"N"],reference)
		version_ch = version_ch.mix(gaps_ch.version)

		exclude_ch = GET_EXCLUDE(meta,reference)
		version_ch = version_ch.mix(exclude_ch.version)

		delly2_ch = DOWNLOAD_DELLY2([:])
		version_ch = version_ch.mix(delly2_ch.version)

	emit:
		executable = delly2_ch.executable
		mappability = map_ch.mappability
		version = version_ch
	}


process GET_MAPPABILITY {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("blacklist.gz"),emit:mappability
	path("version.xml"),emit:version
	path("blacklist.gz.fai")
	path("blacklist.gz.gzi")
script:
	def url19 = "https://gear.embl.de/data/delly/Homo_sapiens.GRCh37.dna.primary_assembly.fa.r101.s501.blacklist.gz"
	def url38 = "https://gear.embl.de/data/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
	def url=(isHg19(reference)?url19:(isHg38(reference)?url38:""))
	
	if(!url.isEmpty())
	"""
	wget -O blacklist.gz "${url}"
	wget -O blacklist.gz.fai "${url}.fai"
	wget -O blacklist.gz.gzi "${url}.gzi"
	
	##########################################################################################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download mappability files</entry>
		<entry key="url">${url}</entry>
	</properties>
	EOF
	"""
 	else
    	"""
 	echo "Undefined build" 1>&2
 	"""
	}



process GET_EXCLUDE {
tag "${file(reference).name}"
afterScript "rm -f jeter.bed jeter2.bed jeter.interval_list"

input:
	val(meta)
	val(reference)
output:
	path("exclude.bed"),emit:bed
	path("version.xml"),emit:version
script:
	def url1=(isHg19(reference)?"https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed":(isHg38(referen,ce)?"https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed":""))
"""
hostname 1>&2
module load ${getModules("jvarkit bedtools picard")}

"""
}



process DOWNLOAD_DELLY2 {
	errorStrategy "retry"
	maxRetries 3
	input:
		val(meta)
	output:
		path("delly"),emit:executable
		path("version.xml"),emit:version
	script:
		def version = getKeyValue(meta,"delly2_version","v1.0.3")
		def url = "https://github.com/dellytools/delly/releases/download/${version}/delly_${version}_linux_x86_64bit"
	"""
	wget -O delly "${url}"
	chmod a+x delly

	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">Download delly</entry>
		<entry key="version">${version}</entry>
		<entry key="url">${url}</entry>
	</properties>
	EOF
	"""
	}


