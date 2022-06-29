include {moduleLoad;getKeyValue} from '../../../modules/utils/functions.nf'
include {VCF_INTER_CASES_CONTROLS_01} from '../../../subworkflows/bcftools/vcf.inter.cases.controls.01.nf'
include {DOWNLOAD_GTF_01} from '../../../modules/gtf/download.gtf.01.nf'
include {BED_CLUSTER_01} from '../../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {SQRT_FILE} from '../../../modules/utils/sqrt.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {WGSELECT_01} from '../../../subworkflows/wgselect/wgselect.01.nf'

workflow {
		BURDEN_CODING(params, params.reference, params.vcf, params.cases, params.controls)
		}

workflow BURDEN_CODING {
	take:
		meta
		reference
		vcf
		cases
		controls
	main:
		version_ch = Channel.empty()
		vcf_inter_ch = VCF_INTER_CASES_CONTROLS_01(meta.plus(["with_tabix":true]),vcf,cases,controls)
		version_ch = version_ch.mix(vcf_inter_ch.version)

		gtf_ch = DOWNLOAD_GTF_01(meta,reference)
		version_ch = version_ch.mix(gtf_ch.version)

		genes_ch = EXTRACT_GENES(meta,reference,gtf_ch.gtf)
		version_ch = version_ch.mix(genes_ch.version)

		copy=genes_ch.genes_bed
		copy.view{"Genes $it"}

                cluster_ch = BED_CLUSTER_01(meta.plus(
                        ["bed_cluster_method":getKeyValue(meta,"bed_cluster_method","--size 1mb")]),
                        reference,
                        genes_ch.genes_bed
                        )
		version_ch = version_ch.mix(cluster_ch.version)

		sqrt_ch = SQRT_FILE(meta.plus(["min_file_split":100]), cluster_ch.output)

		exons_ch = REMOVE_INTRONS(meta,sqrt_ch.clusters.splitText().map{it.trim()}.combine(genes_ch.exons_bed))
		version_ch = version_ch.mix(exons_ch.version)

                file_list_ch = COLLECT_TO_FILE_01([:],exons_ch.bed.splitText().map{it.trim()}.collect())
	
		wgselect_ch = WGSELECT_01(meta,reference,vcf,vcf_inter_ch.cases,vcf_inter_ch.controls,
exons_ch.bed.splitText().map{it.trim()}.map{T->file(T)})

		/*
		RVTESTS_PER_TRANSCRIPT()
		RVTESTS_GROUP()
		PUBLISH()*/
	}

process EXTRACT_GENES {
tag "${gtf.name}"
memory "2g"
input:
	val(meta)
	val(reference)
	path(gtf)
output:
	path("genes.bed"),emit:genes_bed
	path("exons.bed"),emit:exons_bed
	path("version.xml"),emit:version
script:
	def slop=getKeyValue(meta,"genes_slop","50")
"""
hostname 1>&2
${moduleLoad("jvarkit bedtools")}
set -o pipefail

java -jar /${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature,gene_type" "${gtf}" |\
	awk -F '\t' '(\$4=="gene" && \$5=="protein_coding")' |\
	cut -f1,2,3 |\
	bedtools slop -b ${slop} -g "${reference}.fai" |\
	LC_ALL=C sort -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > genes.bed

java -jar /${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature" "${gtf}" |\
	awk -F '\t' '(\$4=="exon")' |\
	cut -f1,2,3 |\
	bedtools slop -b ${slop} -g "${reference}.fai" |\
	LC_ALL=C sort -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > exons.bed

test -s genes.bed
test -s exons.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract protein_coding genes from gtf</entry>
	<entry key="gtf">${gtf}</entry>
	<entry key="slop">${slop}</entry>
	<entry key="bedtools.version">\$(  bedtools --version )</entry>
	<entry key="number.of.genes">\$( wc -l < genes.bed )</entry>
	<entry key="number.of.exons">\$( wc -l < exons.bed )</entry>
</properties>
EOF
"""
}

process REMOVE_INTRONS {
tag "${beds.name}"
input:
	val(meta)
	tuple path(beds),path(exons)
output:
	path("exons.bed.list"),emit:bed
	path("version.xml"),emit:version
"""
hostname 1>&2
set -o pipefail
${moduleLoad("bedtools")}

mkdir BEDS

i=1
cat "${beds}" | while read F
do

	sort -T . -t '\t' -k1,1 -k2,2n "\${F}" |\
		bedtools intersect -u -a "${exons}" -b - > "BEDS/exons.\${i}.bed"
	i=\$((i+1))
done

find \${PWD}/BEDS -type f -name "exons*.bed" > exons.bed.list
test -s exons.bed.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">keep exons from genes bed</entry>
	<entry key="beds">${beds}</entry>
	<entry key="exons">${exons}</entry>
	<entry key="bedtools.version">\$(  bedtools --version )</entry>
</properties>
EOF
"""
}

