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

include {SAMTOOLS_SAMPLES02} from '../../modules/samtools/samtools.samples.02.nf'
include {moduleLoad;assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {CONCAT_FILES_01} from '../../modules/utils/concat.files.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SAMTOOLS_DEPTH_VCF_01 as ST_DEPTH} from '../../modules/samtools/samtools.depth.vcf.01.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'


workflow SAMTOOLS_DEPTH_VCF_01 {
	take:
		meta
		reference
		references
		bams
		bed
	main:
		version_ch = Channel.empty()
		
		samples_ch = SAMTOOLS_SAMPLES02(
			meta.plus(["with_header":"true"]),
			["reference":reference,"references":references,"bams":bams]
			)
		version_ch = version_ch.mix(samples_ch.version)

		dp_ch = ST_DEPTH(meta,samples_ch.output.splitCsv(header:true,sep:"\t").
			map{T->T.plus("bed":bed)}
			)
		version_ch = version_ch.mix(dp_ch.version)


		to_file_ch = COLLECT_TO_FILE_01(meta, dp_ch.output.map{T->T[1]}.collect())
		version_ch= version_ch.mix(to_file_ch.version)

		sqrt_ch = SQRT_FILE(meta.plus(["min_file_split":100]), to_file_ch.output)
		version_ch= version_ch.mix(sqrt_ch.version)

		each_cluster = sqrt_ch.output.splitText().map{T->file(T.trim())}

		merge0_ch = MERGE_BCFS0(meta,each_cluster)
		version_ch = version_ch.mix(merge0_ch.version)

		merge1_ch = MERGE_BCFS1(meta,reference,merge0_ch.vcf.collect())
		version_ch = version_ch.mix(merge1_ch.version)
		
		version_ch = MERGE_VERSION(meta, "Samtools depth", "samtools depth", version_ch.collect())
	emit:
		version = version_ch
		vcf = merge1_ch.vcf
		index = merge1_ch.index
	}


process MERGE_BCFS0 {
tag "${vcfs.name}"
cpus 10
input:
	val(meta)
	path(vcfs)
output:
	path("merged.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

bcftools merge --threads ${task.cpus} --file-list "${vcfs}" --info-rules DP:sum,MINDP:min,MAXDP:max --no-version -O b -o merged.bcf
bcftools index merged.bcf



##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">bcftools merge</entry>
        <entry key="vcfs">${vcfs}</entry>
        <entry key="bcftools.version">\$(bcftools version | head -n2 | paste -s )</entry>
</properties>
EOF
"""
}

process MERGE_BCFS1 {
tag "N=${L.size()}"
cpus 10
memory "7g"
afterScript "rm -f jeter.bcf jeter2.bcf jeter.list"
input:
	val(meta)
	val(reference)
	val(L)
output:
	path("${params.prefix?:""}depth.base.bcf"),emit:vcf
	path("${params.prefix?:""}depth.base.bcf.csi"),emit:index
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail

cat << EOF > jeter.list
${L.join("\n")}
EOF

bcftools merge --threads ${task.cpus} \
	--info-rules DP:sum,MINDP:min,MAXDP:max \
	--no-version \
	-O u \
	--file-list jeter.list \
	-o jeter.bcf

bcftools sort --max-mem "${task.memory.giga}G" -T . -O b -o jeter2.bcf jeter.bcf
mv jeter2.bcf jeter.bcf

bcftools norm --check-ref s --fasta-ref "${reference}" -O b -o "${params.prefix?:""}depth.base.bcf" jeter.bcf

bcftools index "${params.prefix?:""}depth.base.bcf"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">bcftools merge</entry>
        <entry key="vcfs.count">${L.size()}</entry>
        <entry key="bcftools.version">\$(bcftools version | head -n2 | paste -s )</entry>
</properties>
EOF
"""
}

