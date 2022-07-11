include {moduleLoad;assertFileExists;isBlank} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {DOWNLOAD_SOMALIER} from '../../modules/somalier/somalier.download.nf'
include {SOMALIER_DOWNLOAD_SITES} from '../../modules/somalier/somalier.download.sites.nf'
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'


workflow SOMALIER_BAMS_01 {
	take:
		meta
		reference
		bams
		pedigree
	main:
		assertFileExists(reference,"reference must be defined")
		version_ch = Channel.empty()

		exe_ch = DOWNLOAD_SOMALIER(meta)
		version_ch = version_ch.mix(exe_ch.version)

		sites_ch = SOMALIER_DOWNLOAD_SITES(meta,reference)
		version_ch = version_ch.mix(sites_ch.version)
		
		ch1 = SAMTOOLS_SAMPLES01([:],reference,bams)
		version_ch = version_ch.mix(ch1.version)

		ch2 = ch1.output.splitCsv(header:false,sep:'\t')		

		ch3 = EXTRACT_BAM(meta, reference, exe_ch.executable , sites_ch.vcf , ch2)
		version_ch = version_ch.mix(ch3.version)

		somalier_ch = RELATE_SOMALIER(meta,reference, exe_ch.executable,ch3.output.collect(), pedigree)
		version_ch = version_ch.mix(somalier_ch.version)

		version_ch = MERGE_VERSION(meta, "somalier", "somalier on bams", version_ch.collect())
	emit:
		version = version_ch
		zip = somalier_ch.zip
	}


process EXTRACT_BAM {
	tag "${sample}"
	memory '2g'

	input:
		val(meta)
		val(reference)
		val(somalier)
		val(sites)
		tuple val(sample),val(bam)
	output:
		path("extracted/${sample}.somalier"),emit:output
		path("version.xml"),emit:version
	script:
	"""
	hostname 1>&2
	mkdir extracted
	${somalier} extract -d extracted --sites "${sites}" -f "${reference}" "${bam}"
	
	test -s "extracted/${sample}.somalier"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">somalier extract</entry>
        <entry key="sample">${sample}</entry>
        <entry key="bam">${bam}</entry>
</properties>
EOF
	"""
	}


process RELATE_SOMALIER {
tag "N=${L.size()}"
afterScript "rm -rf extracted TMP"
input:
	val(meta)
	val(reference)
	val(somalier)
	val(L)
	path(pedigree)
output:
	path("${meta.prefix?:""}somalier.bams.zip"),emit:zip
	path("version.xml"),emit:version
script:
	def prefix = meta.prefix?:""

"""
hostname 1>&2
set -o pipefail
set -x

mkdir TMP

mkdir "${prefix}somalier.bams"

${somalier} relate --output-prefix=${prefix}somalier.bams/${prefix}bams \
	${pedigree.name.equals("NO_FILE")?"":"-p '${pedigree}'"} \
	${L.join(" ")}

# may not exist
touch "${prefix}somalier.bams/${prefix}bams.groups.tsv"
zip -9 -r "${prefix}somalier.bams.zip" "${prefix}somalier.bams"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">run somalier on bams</entry>
        <entry key="bams.count">${L.size()}</entry>
        <entry key="pedigree">${pedigree}</entry>
</properties>
EOF
"""
}
