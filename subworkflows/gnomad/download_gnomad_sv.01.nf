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

include {moduleLoad;getVersionCmd;isHg19;isHg38} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

def gazoduc = gazoduc.Gazoduc.getInstance(params)

gazoduc.build("gnomad_sv_grch37_bed_url","https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz").
	menu("Gnomad").
        desc("URL for gnomad SV on build grch37.").
        put()


workflow DOWNLOAD_GNOMAD_SV_01 {
     take:
        meta /* meta */
        reference /* indexed fasta reference either hg38 of hg39*/
     main:
        version_ch = Channel.empty()

	ch0 = DOWNLOAD_BED(meta)
	version_ch = version_ch.mix(ch0.version)

	if(isHg19(reference)) {
		ch1 = RENAME_CONTIGS_HG19(meta, ch0.output, reference)
		version_ch = version_ch.mix(ch1.version)

		dest_bed = ch1.bed
		dest_idx = ch1.index
		}
	else if(isHg38(reference)) {
		ch2 = LIFT_GNOMAD_SV_TO_HG38(meta, ch0.output, reference)
		version_ch = version_ch.mix(ch2.version)

		dest_bed = ch2.bed
		dest_idx = ch2.index
		}
	else
		{
		ch3 = OTHER_REF(meta, ch0.output)

		version_ch = version_ch.mix(ch3.version)

		dest_bed = ch3.bed
		dest_idx = ch3.index
		}
		
	version_ch = MERGE_VERSION(meta, "GnomadSV", "GnomadSV",version_ch.collect())

	emit:
		output = dest_bed
		index = dest_idx
		version = version_ch
	}

process DOWNLOAD_BED {
tag "${meta.gnomad_sv_bed_url}"
input:
	val(meta)
output:
	path("gnomad.sv.37.bed.gz"),emit:output
	path("version.xml"),emit:version
script:
	def url = meta.gnomad_sv_grch37_bed_url
"""
wget -O "gnomad.sv.37.bed.gz" "${url}" 

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">Download gnomad sv</entry>
        <entry key="gnomad.sv.url"><a>${url}</a></entry>
        <entry key="versions">${getVersionCmd("wget")}</entry>
</properties>
EOF
"""
}

process RENAME_CONTIGS_HG19 {
tag "${meta.gnomad_sv_bed_url}"
input:
	val(meta)
	path(bed)
	val(reference)
output:
	path("gnomad.37.sv.bed.gz"),emit:bed
	path("gnomad.37.sv.bed.gz.tbi"),emit:index
	path("version.xml"),emit:version
script:
"""
${moduleLoad("jvarkit htslib")}
set -o pipefail
mkdir -p TMP

gunzip -c "${bed}" | head -n1  > TMP/gnomad.37.sv.bed


gunzip -c "${bed}" |\
	tail -n+2 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n >> TMP/gnomad.37.sv.bed

bgzip -f TMP/gnomad.37.sv.bed
tabix --comment '#'  -p bed -f TMP/gnomad.37.sv.bed.gz

mv TMP/gnomad.37.sv.bed.gz ./
mv TMP/gnomad.37.sv.bed.gz.tbi ./

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">fix chromosomes for hg19</entry>
</properties>
EOF
"""
}


process LIFT_GNOMAD_SV_TO_HG38 {
tag "${meta.gnomad_sv_bed_url}"
input:
	val(meta)
	path(bed)
	val(reference)
output:
	path("gnomad.38.sv.bed.gz"),emit:bed
	path("gnomad.38.sv.bed.gz.tbi"),emit:index
	path("version.xml"),emit:version
script:
	def url1 = "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19Tohg38.over.chain.gz"
"""
${moduleLoad("ucsc jvarkit htslib")}
set -o pipefail
mkdir -p TMP

wget -O TMP/hg19Tohg38.over.chain.gz "${url1}"
gunzip -f  TMP/hg19Tohg38.over.chain.gz


gunzip -c "${bed}" | head -n1  > TMP/gnomad.38.sv.bed


gunzip -c "${bed}" |\
	tail -n+2 |\
	sed 's/^/chr/' |\
	liftover -bedPlus=3 stdin TMP/hg19Tohg38.over.chain stdout TMP/unmapped.bed |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n >> TMP/gnomad.38.sv.bed

bgzip -f TMP/gnomad.38.sv.bed
tabix --comment '#'  -p bed -f TMP/gnomad.38.sv.bed.gz

mv TMP/gnomad.38.sv.bed.gz ./
mv TMP/gnomad.38.sv.bed.gz.tbi ./

#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">Download gnomad sv and liftover 38</entry>
</properties>
EOF
"""
}

process OTHER_REF {
input:
	val(meta)
	path(bed)
output:
	path("gnomad.x.sv.bed.gz"),emit:bed
	path("gnomad.x.sv.bed.gz.tbi"),emit:index
	path("version.xml"),emit:version
script:
	def url2 = meta.gnomad_sv_grch37_bed_url
"""
${moduleLoad("htslib")}
gunzip -c "${bed}" | head -n 1 > gnomad.x.sv.bed
bgzip -f gnomad.x.sv.bed
tabix --comment '#'  -p bed -f  gnomad.x.sv.bed.gz


#########################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="Name">${task.process}</entry>
        <entry key="description">Download gnomad sv but undefined ref</entry>
</properties>
EOF
"""
}
