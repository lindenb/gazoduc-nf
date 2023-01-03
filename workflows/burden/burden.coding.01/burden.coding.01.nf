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
nextflow.enable.dsl=2

include {getVersionCmd;runOnComplete;moduleLoad;getKeyValue;hasFeature} from '../../../modules/utils/functions.nf'
include {VALIDATE_CASE_CONTROL_PED_01} from '../../../modules/pedigree/validate.case.ctrl.pedigree.01.nf'
include {VCF_INTER_CASES_CONTROLS_01} from '../../../subworkflows/bcftools/vcf.inter.cases.controls.01.nf'
include {DOWNLOAD_GTF_01} from '../../../modules/gtf/download.gtf.01.nf'
include {BED_CLUSTER_01} from '../../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {SQRT_FILE} from '../../../modules/utils/sqrt.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
include {WGSELECT_01} from '../../../subworkflows/wgselect/wgselect.01.nf'
include {PEDIGREE_FOR_RVTESTS} from '../../../modules/rvtests/rvtests.cases.controls.ped.01.nf'
include {RVTESTS01_VCF_01} from '../../../modules/rvtests/rvtests.vcf.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {PIHAT_CASES_CONTROLS_01} from '../../../subworkflows/pihat/pihat.cases.controls.01.nf'

params.reference=""
params.pedigree=""
params.vcf=""
params.disableFeatures="";
params.help=false
params.bed= "NO_FILE"
params.annot_method="vep"

if(params.help) {
  log.info"""
## About

Burden for coding regions.

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
  * --bed <file> optional bed file to limit the analysis to the genes overlapping a  bed file.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
        --publishDir output \\
        --prefix "analysis." \\
        --reference /path/to/reference.fasta \\
        --vcf /path/to/my.vcf.gz \\
        --pedigree /path/to/input.ped
```

## Workflow

![workflow](./workflow.svg)

"""
exit 0
}

workflow {
		BURDEN_CODING(params, params.reference, params.vcf, file(params.pedigree),file(params.bed))
		}

workflow BURDEN_CODING {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:
		to_zip = Channel.empty()
		version_ch = Channel.empty()
		
		ped_ch = VALIDATE_CASE_CONTROL_PED_01(meta,pedigree)
		version_ch = version_ch.mix(ped_ch.version)

		vcf_inter_ch = VCF_INTER_CASES_CONTROLS_01(meta.plus(["with_tabix":true]),vcf,
					ped_ch.cases_list, ped_ch.controls_list )
		version_ch = version_ch.mix(vcf_inter_ch.version)

		if(hasFeature(meta,"pihat")) {
			pihat = PIHAT_CASES_CONTROLS_01(meta,reference,file(vcf),ped_ch.cases_list,ped_ch.controls_list)
			version_ch = version_ch.mix(pihat.version)
			to_zip = to_zip.mix(pihat.pihat_png)
			to_zip = to_zip.mix(pihat.pihat_pdf)
			to_zip = to_zip.mix(pihat.removed_samples)
			to_zip = to_zip.mix(pihat.plink_genome)

			rebuild_ch = REBUILD_PEDIGREE(meta,pedigree,pihat.cases,pihat.controls)
			version_ch = version_ch.mix(rebuild_ch.version)

			new_ped_ch = rebuild_ch.pedigree
			}
		else
			{
			rebuild_ch = REBUILD_PEDIGREE(meta,pedigree,ped_ch.cases_list,ped_ch.controls_list)
			version_ch = version_ch.mix(rebuild_ch.version)

			new_ped_ch = rebuild_ch.pedigree
			}


		gtf_ch = DOWNLOAD_GTF_01(meta,reference)
		version_ch = version_ch.mix(gtf_ch.version)

		genes_ch = EXTRACT_GENES(meta,reference,gtf_ch.gtf,bed)
		version_ch = version_ch.mix(genes_ch.version)

                cluster_ch = BED_CLUSTER_01(meta.plus(
                        ["bed_cluster_method":getKeyValue(meta,"bed_cluster_method","--size 5mb")]),
                        reference,
                        genes_ch.genes_bed
                        )
		version_ch = version_ch.mix(cluster_ch.version)

		sqrt_ch = SQRT_FILE(meta.plus(["min_file_split":100]), cluster_ch.output)

		exons_ch = REMOVE_INTRONS(meta,sqrt_ch.clusters.splitText().map{it.trim()}.combine(genes_ch.exons_bed))
		version_ch = version_ch.mix(exons_ch.version)

                file_list_ch = COLLECT_TO_FILE_01([:],exons_ch.bed.splitText().map{it.trim()}.collect())
	
		wgselect_ch = WGSELECT_01(meta,reference,vcf,new_ped_ch,exons_ch.bed.splitText().map{it.trim()}.map{T->file(T)})
		version_ch = version_ch.mix(wgselect_ch.version.first())
		to_zip = to_zip.mix(wgselect_ch.variants_list)

		rvtests_ped_ch = PEDIGREE_FOR_RVTESTS(meta,new_ped_ch)
		version_ch = version_ch.mix(rvtests_ped_ch.version)

		assoc_ch = RVTESTS01_VCF_01(meta, reference,wgselect_ch.vcfs.splitText().map{it.trim()}, rvtests_ped_ch.pedigree )
		version_ch = version_ch.mix(assoc_ch.version.first())

		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference,file(vcf).name,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden coding", "Burden coding ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
		PUBLISH(to_zip.collect())
	}

process EXTRACT_GENES {
tag "${gtf.name}"
memory "2g"
input:
	val(meta)
	val(reference)
	path(gtf)
	path(bed)
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

java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature,gene_type" "${gtf}" |\
	awk -F '\t' '(\$4=="gene" && \$5=="protein_coding")' |\
	cut -f1,2,3 |\
	bedtools slop -b ${slop} -g "${reference}.fai" |\
	LC_ALL=C sort -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > genes.bed

java -jar \${JVARKIT_DIST}/gtf2bed.jar -R "${reference}" --columns "gtf.feature" "${gtf}" |\
	awk -F '\t' '(\$4=="exon")' |\
	cut -f1,2,3 |\
	bedtools slop -b ${slop} -g "${reference}.fai" |\
	LC_ALL=C sort -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > exons.bed


if [ ! -z "${bed.name.equals("NO_FILE")?"":"Y"}" ] ; then

	LC_ALL=C sort -T . -t '\t' -k1,1 -k2,2n "${bed}" | cut -f1,2,3 > jeter.bed
	test -s jeter.bed

	bedtools intersect -u -a genes.bed -b jeter.bed |\
	LC_ALL=C sort -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n > jeter2.bed
	mv jeter2.bed genes.bed

	bedtools intersect -u -a exons.bed -b jeter.bed |\
	LC_ALL=C sort -S ${task.memory.kilo} -T . -t '\t' -k1,1 -k2,2n > jeter2.bed
	mv jeter2.bed exons.bed

fi

test -s genes.bed
test -s exons.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">extract protein_coding genes from gtf</entry>
	<entry key="gtf">${gtf}</entry>
	<entry key="slop">${slop}</entry>
	<entry key="versions">${getVersionCmd("bedtools awk jvarkit/gtf2bed")}</entry>
	<entry key="number.of.genes">\$( wc -l < genes.bed )</entry>
	<entry key="number.of.exons">\$( wc -l < exons.bed )</entry>
	<entry key="bed">${bed}</entry>
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
	<entry key="bedtools.version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}

process REBUILD_PEDIGREE {
executor "local"
input:
	val(meta)
	path(pedigree)
	path(cases)
	path(controls)
output:
	path("case_control.ped"),emit:pedigree
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail

sort -T . -t '\t' -k2,2 "${pedigree}" > tmp1.ped
cat  "${cases}" "${controls}" | sort | uniq > tmp2.ped

join -t '\t' -1 2 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6' tmp1.ped  tmp2.ped |\
	sort -T . -t '\t' -k1,1 -k2,2  > case_control.ped
test -s case_control.ped
rm tmp1.ped tmp2.ped



###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">make pedigree after samples have been removed.</entry>
	<entry key="pedigree">${pedigree.toRealPath()}</entry>
	<entry key="cases">${cases.toRealPath()}</entry>
	<entry key="controls">${controls.toRealPath()}</entry>
	<entry key="samples.removed">\$(comm -23 <(cut -f 2 "${pedigree}" | sort)   <(cut -f 2 "case_control.ped" | sort)  | paste -s -d ' ' )</entry>
</properties>
EOF
"""
}

process PUBLISH {
tag "N=${files.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
        val(files)
output:
        path("${params.prefix?:""}archive.zip")
when:
        !params.getOrDefault("publishDir","").trim().isEmpty()
script:
        prefix = params.getOrDefault("prefix","")
"""

mkdir "${prefix}archive"

cat << EOF | while read F ; do ln -s "\${F}" "./${prefix}archive/" ; done
${files.join("\n")}
EOF

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}

runOnComplete(workflow);
