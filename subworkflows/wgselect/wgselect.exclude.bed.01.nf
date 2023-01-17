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
include {isHg19;isHg38;moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include { SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'

def gazoduc = gazoduc.Gazoduc.getInstance(params)

gazoduc.
    make("wgselect_with_rmsk",true).
    menu("wgselect").
    description("remove variant overlapping repeat masker UCSC regions").
    setBoolean().
    put()

gazoduc.
    make("wgselect_with_encodeExclude",true).
    menu("wgselect").
    description("remove variant overlapping excluded ENCODe regions").
    setBoolean().
    put()

gazoduc.
    make("wgselect_with_lcr",true).
    menu("wgselect").
    description("remove variant overlapping low complexity regions").
    setBoolean().
    put()

gazoduc.
    make("wgselect_with_simpleRepeats",true).
    menu("wgselect").
    description("remove variant overlapping simple repeats regions").
    setBoolean().
    put()

workflow WGSELECT_EXCLUDE_BED_01 {
	take:
		meta
		reference
	main:
		version_ch = Channel.empty()
		to_merge_ch = Channel.empty()

		gaps_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"N","MAX_TO_MERGE":"1"],reference)
		to_merge_ch = to_merge_ch.mix(gaps_ch.bed)

		if(meta.wgselect_with_rmsk as boolean) {
			rmsk_ch = RMSK(meta,reference)
			version_ch = version_ch.mix(rmsk_ch.version)
			to_merge_ch = to_merge_ch.mix(rmsk_ch.bed)
			}

		if(meta.wgselect_with_encodeExclude as boolean) {
			x2_ch = EXCLUDE_ENCODE(meta,reference)
			version_ch = version_ch.mix(x2_ch.version)
			to_merge_ch = to_merge_ch.mix(x2_ch.bed)
			}


		if(meta.wgselect_with_lcr as boolean) {
			x3_ch = LOW_COMPLEXITY_REGIONS(meta,reference)
			version_ch = version_ch.mix(x3_ch.version)
			to_merge_ch = to_merge_ch.mix(x3_ch.bed)
			}

		if(meta.wgselect_with_simpleRepeats as boolean) {
			x4_ch = SIMPLE_REPEATS(meta,reference)
			version_ch = version_ch.mix(x4_ch.version)
			to_merge_ch = to_merge_ch.mix(x4_ch.bed)
			}
		all_x_ch = MERGE_REGIONS(meta,to_merge_ch.collect())
		version_ch = version_ch.mix(all_x_ch.version)

		version_ch = MERGE_VERSION(meta,"blacklisted","blacklisted regions",version_ch.collect())
	emit:
		version = version_ch.version
		bed = all_x_ch.bed
	}




process RMSK {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("rmsk.bed"), emit:bed
	path("version.xml"), emit:version
script:
	def base="https://hgdownload.cse.ucsc.edu/goldenPath/_DB_/database/rmsk.txt.gz"
	def url =	(isHg19(reference)?base.replaceAll("_DB_","hg19"):
			(isHg38(reference)?base.replaceAll("_DB_","hg38"):
			""));
	if(url.isEmpty()) throw new IllegalArgumentException("undefined rmsk.url for ${reference}");

"""
hostname 1>&2
set -o pipefail
${moduleLoad("jvarkit")}

wget -O - "${url}" |\
	gunzip -c |\
	cut -f6-8 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
		LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > rmsk.bed

test -s rmsk.bed

#########################################################"
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">load bed from UCSC : repeat masked regions. To disable add 'rmsk' to '<code>--disableFeatures</code>'</entry>
	<entry key="url"><a>${url}</a></entry>
	<entry key="jvarkit.bedrename.version">\$( java -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
</properties>
EOF
"""
}

process EXCLUDE_ENCODE {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("excude.encode.bed"),emit:bed
	path("version.xml"),emit:version
script:
	// https://www.biostars.org/p/171354/
	def url =	(isHg19(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz":
			(isHg38(reference)?"https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz?raw=true":
			""));
	if(url.isEmpty()) throw new IllegalArgumentException("undefined encode blakclist for ${reference}");

"""
hostname 1>&2
set -o pipefail
${moduleLoad("jvarkit")}

wget -O - "${url}" |\
	gunzip -c |\
	cut -f1-3 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > excude.encode.bed

test -s excude.encode.bed


#########################################################"
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">load excludedRegions from encode. To disable add 'encodeExclude' to '<code>--disableFeatures</code>'</entry>
	<entry key="url"><a>${url}</a></entry>
	<entry key="jvarkit.bedrename.version">\$( java -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
</properties>
EOF
"""
}


process LOW_COMPLEXITY_REGIONS {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("lcr.bed"), emit:bed
	path("version.xml"), emit:version
script:
	def url =	(isHg19(reference)?"https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs37d5.bed.gz?raw=true":
			(isHg38(reference)?"https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true":
			""));
	if(url.isEmpty()) throw new IllegalArgumentException("undefined encode blakclist for ${reference}");

"""
hostname 1>&2
set -o pipefail
${moduleLoad("jvarkit")}
		
wget -O - "${url}" |\
	gunzip -c |\
	cut -f1-3 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n > lcr.bed

test -s lcr.bed

#########################################################"
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">download lowComplexity Regions. To disable add 'lcr' to '<code>--disableFeatures</code>'</entry>
	<entry key="url"><a>${url}</a></entry>
	<entry key="jvarkit.bedrename.version">\$( java -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
</properties>
EOF
"""
}

process SIMPLE_REPEATS {
tag "${file(reference).name}"
input:
	val(meta)
	val(reference)
output:
	path("simple_repeats.bed"), emit:bed
	path("version.xml"), emit:version
script:
	def url =	(isHg19(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz":
			(isHg38(reference)?"https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz":
			""));

	if(url.isEmpty()) throw new IllegalArgumentException("undefined encode blakclist for ${reference}");
"""
hostname 1>&2
set -o pipefail
${moduleLoad("jvarkit")}

wget -O - "${url}" |\
	gunzip -c |\
	cut -f2-4 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
		LC_ALL=C sort -t '\t' -T . -k1,1 -k2,2n > simple_repeats.bed


test -s simple_repeats.bed

#########################################################"
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">download simple repeats. To disable add 'simpleRepeats' to '<code>--disableFeatures</code>'</entry>
	<entry key="url"><a>${url}</a></entry>
	<entry key="jvarkit.bedrename.version">\$( java -jar \${JVARKIT_DIST}/bedrenamechr.jar --version)</entry>
</properties>
EOF
"""
}


process MERGE_REGIONS {
tag "N=${L.size()}"
input:
	val(meta)
	val(L)
output:
	path("bad_regions.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools")}
set -o pipefail
LC_ALL=C sort --merge -T . -k1,1 -k2,2n ${L.join(" ")} |\
	cut -f1-3 |\
	bedtools merge > bad_regions.bed


#########################################################"
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge all ${L.size()} blacklisted bed(s)</entry>
	<entry key="bedtools.version">\$(bedtools --version)</entry>
</properties>
EOF
"""
}


