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

include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;isHg19;isHg38;runOnComplete;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VCF_TO_BED} from '../../modules/bcftools/vcf2bed.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {SAMPLES_IN_VCF_01} from '../../modules/bcftools/samples.in.vcf.01.nf'


workflow VCF_DUPHOLD_01 {
	take:
		meta
		reference
		bams
		vcf
		snps /** small indels , snps */
	main:
		version_ch  = Channel.empty()

		bams_ch = SAMTOOLS_SAMPLES_01(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false],reference,file("NO_FILE"),bams)
		version_ch = version_ch.mix(bams_ch.version)
	
		vcf2bed_ch = VCF_TO_BED(meta, vcf)
		version_ch = version_ch.mix(vcf2bed_ch.version)

		contig_bed_ch = vcf2bed_ch.bed.splitCsv(header:false,sep:'\t').map{T->[T[0],T[3]]}

		snvcf_ch = SAMPLES_IN_VCF_01(meta, vcf)
		version_ch = version_ch.mix(snvcf_ch.version)

		to_concat_ch = Channel.empty()

		dispatch_ch = DISPATCH_VCF(meta, contig_bed_ch)
		version_ch = version_ch.mix(dispatch_ch.version)
		to_concat_ch = to_concat_ch.mix(dispatch_ch.remain)

		join_ch = JOIN_BAM_VCFS(meta, bams_ch.output, snvcf_ch.samples)
		version_ch = version_ch.mix(join_ch.version)

		exec_ch = DOWNLOAD_DUPHOLD(meta)
		version_ch = version_ch.mix(exec_ch.version)
		
		call_ch= APPLY_DUPHOLD(meta, reference, exec_ch.executable, snps, join_ch.bams, dispatch_ch.output.splitText().map{T->T.trim()} )
		version_ch = version_ch.mix(call_ch.version)
		to_concat_ch = to_concat_ch.mix(call_ch.vcf)

		x3_ch = COLLECT_TO_FILE_01([:],to_concat_ch.collect())
		version_ch = version_ch.mix(x3_ch.version)


		concat_ch = BCFTOOLS_CONCAT_01([:], x3_ch.output)
		version_ch = version_ch.mix(concat_ch.version)
		
		
		version_ch = MERGE_VERSION(meta, "Duphold", "Duphold",version_ch.collect())		

	emit:
		vcf = concat_ch.vcf
		index = concat_ch.index
		version= version_ch
	}


process DISPATCH_VCF {
tag "${contig} ${file(vcf).name}"
input:
	val(meta)
	tuple val(contig),val(vcf)
output:
	path("vcfs.list"),emit:output
	path("remain.${contig}.bcf"),emit:remain
	path("version.xml"),emit:version
script:
	def method = meta.duphold_split_method?:"--variants-count 100"
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}
set -o pipefail
mkdir -p OUT

bcftools view -i 'SVTYPE=="DEL" || SVTYPE=="DUP"' "${vcf}" "${contig}" |\
	java -jar \${JVARKIT_DIST}/vcfsplitnvariants.jar ${method} -o \${PWD}/OUT/split.${contig}.

find OUT/ -type f -name "*.vcf.gz" | while read F
do
	bcftools index -t --force "\${F}"
done

find \${PWD}/OUT/ -type f -name "*.vcf.gz" > vcfs.list

# remaining
bcftools view -e 'SVTYPE=="DEL" || SVTYPE=="DUP"' "${vcf}" "${contig}" -O b -o "remain.${contig}.bcf"
bcftools index "remain.${contig}.bcf"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">split VCF per contig and using jvarkit/vcfsplitnvariants</entry>
        <entry key="method">${method}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="contig">${contig}</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/vcfsplitnvariants")}</entry>
</properties>
EOF
"""
}

process JOIN_BAM_VCFS {
executor "local"
tag "${samples} ${bams}"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(bams)
	path(samples)
output:
	path("join.bams.list"),emit:bams
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
mkdir -p TMP

cut -f1,3 "${bams}" | sort -t '\t' -k1,1 -T TMP > TMP/jeter.a

sort -t '\t' -k1,1 -T TMP "${samples}"> TMP/jeter.b

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter.a TMP/jeter.b > join.bams.list

test -s join.bams.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join vcf samples and bams</entry>
        <entry key="samples">${samples}</entry>
        <entry key="bams">${bams}</entry>
</properties>
EOF
"""
}


process DOWNLOAD_DUPHOLD {
input:
	val(meta)
output:
	path("duphold"),emit:executable
	path("version.xml"),emit:version
script:
	def version = meta.duphold_version?:"0.2.3"
	def url = "https://github.com/brentp/duphold/releases/download/v${version}/duphold"

"""
hostname 1>&2
wget -O duphold "${url}"
chmod +x duphold

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download duphold</entry>
        <entry key="version">${version}</entry>
        <entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

/*
process DOWNLOAD_GNOMAD {
input:
	val(meta)
	val(reference)
output:
	path("gnomad.sites.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz":(isHg38(reference)?"":"")
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools")}
set -o pipefail


wget -O - "${url}" |\
	bcftools view |\
	java -jar \${JVARKIT_DIST}/vcfsetdict.jar -R "${reference}"  --onNotFound SKIP |\
	bcftools sort -T . -O b -o gnomad.sites.bcf

bcftools index gnomad.sites.bcf

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download gnomad</entry>
        <entry key="versions">${getVersionCmd("bcftools jvarkit/vcfsetdict")}</entry>
</properties>
EOF
"""
} */



process APPLY_DUPHOLD {
tag "${vcf.name}"
afterScript "rm -rf TMP"
memory "10g"
cpus 4
input:
	val(meta)
	val(reference)
	val(duphold)
	path(snps)//TODO use it
	path(bams)
	path(vcf)
output:
	path("duphold.bcf"),emit:vcf
	path("version.xml"),emit:version
script:

"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir TMP

bcftools view --threads "${task.cpus}" -O b -o TMP/jeter.bcf "${vcf}"

cat "${bams}" | while read B
do

echo "#Processing \${B}" 1>&2

${duphold} --vcf "TMP/jeter.bcf" \
	--threads "${task.cpus}"  \
	--bam "\${B}" \
	--fasta "${reference}" \
	--output TMP/jeter2.bcf 1>&2

mv -v TMP/jeter2.bcf  TMP/jeter.bcf

done

bcftools view --threads "${task.cpus}" -O b -o "duphold.bcf" TMP/jeter.bcf
bcftools index --threads "${task.cpus}" "duphold.bcf"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">apply duphold</entry>
        <entry key="reference">${reference}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="bams">${bams}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}


process ZIP {
executor "local"
input:
	val(meta)
	path(concordances)
	path(version)
	path(html)
	path(pdf)
output:
	path("${meta.prefix?:""}archive.zip"),emit:zip
script:
"""
cp "${concordances}" jeter.list
echo "${version}" >> jeter.list
echo "${html}" >> jeter.list
echo "${html}" >> jeter.list
echo "${pdf}" >> jeter.list

zip -9 -@ -j "${meta.prefix?:""}archive.zip" < jeter.list
rm jeter.list
"""
}
