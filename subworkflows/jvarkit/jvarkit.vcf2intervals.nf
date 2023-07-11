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
log.info("parsing subworkflows/jvarkit/jvarkit.vcf2intervals.nf ....")

include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {JVARKIT_VCF_TO_INTERVALS_01 as VCF_TO_INTERVALS} from '../../modules/jvarkit/jvarkit.vcf2intervals.01.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'



/**
 * call jvarkit vcf2intervals for one or more vcf
 *
 */
workflow JVARKIT_VCF_TO_INTERVALS_01 {
	take:
		meta
		vcf /* path to vcf, or file with .list suffix */
		bed /* limit to that BED or NO_FILE */
	main:
		if(!meta.containsKey("with_header")) throw new RuntimeException("meta.width_header missing")
		if(!meta.containsKey("distance")) throw new RuntimeException("meta.distance missing")
		if(!meta.containsKey("min_distance")) throw new RuntimeException("meta.min_distance missing")
	
		version_ch = Channel.empty()
		bed_ch = VCF_TO_BED([with_header:false],vcf)
		version_ch = version_ch.mix(bed_ch.version)
	
		bed2ch = INTERSECT_BED(bed_ch.bed, bed)
		version_ch = version_ch.mix(bed2ch.version)

		interval_vcf = bed2ch.bed.splitCsv(header:false,sep:'\t').
			map{T->[vcf:file(T[3]),interval:T[0]+":"+((T[1] as int)+1)+"-"+T[2]]}

		rch = VCF_TO_INTERVALS(meta,interval_vcf)
		version_ch = version_ch.mix(rch.version)

		concat_ch = CONCAT_FILES_01([suffix:".bed", concat_n_files:50,downstream_cmd:""] , rch.bed.collect())
		version_ch = version_ch.mix(concat_ch.version)

		version_ch = MERGE_VERSION("jvarkit2intervals",version_ch.collect())

	emit:
		bed = concat_ch.output
		version = version_ch
	}

process INTERSECT_BED {
tag "${vcfbed.name} / ${userbed.name}"
executor "local"
afterScript "rm -rf TMP"
input:
	path(vcfbed)
	path(userbed)
output:
	path("intersect.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bedtools")}
set -o pipefail

mkdir TMP

if ${userbed.name.equals("NO_FILE")} ; then
	cp "${vcfbed}" intersect.bed
else

	sort -t '\t' -T TMP -k1,1 -k2,2n "${vcfbed}" >  TMP/a.bed
	sort -t '\t' -T TMP -k1,1 -k2,2n "${userbed}" | cut -f1,2,3 | bedtools merge >  TMP/b.bed
	bedtools intersect -a TMP/a.bed -b TMP/b.bed |\
		sort -t '\t' -T TMP -k1,1 -k2,2n > intersect.bed
fi

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">intersection between vcf bed and user bed</entry>
	<entry key="version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}


