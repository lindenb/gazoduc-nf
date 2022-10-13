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

include {isHg19;runOnComplete;moduleLoad;getKeyValue;hasFeature;getVersionCmd} from '../../../modules/utils/functions.nf'
include {BURDEN_SAMPLE_WGSELECT_PART_01}  from '../../../subworkflows/burden/burden.samples.wgselect.part.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {RVTESTS_REHEADER_01} from '../../../modules/rvtests/rvtests.reheader.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'

params.reference=""
params.pedigree=""
params.vcf=""
params.help=false
params.bed_cluster_method = " --size 1mb "
params.bed=""

if(params.help) {
  log.info"""
## About

Burden for 1st intron.

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
  * --bed <file> CHROM/START/END/GENE_NAME.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume worfklow.nf \\
        --publishDir output \\
        --prefix "analysis." \\
        --reference /path/to/reference.fasta \\
        --vcf /path/to/my.vcf.gz \\
        --pedigree /path/to/input.ped \
```

## Workflow

![workflow](./workflow.svg)

"""
exit 0
}

workflow {
		burden_ch = BURDEN_PAIRS(params, params.reference, params.vcf, file(params.pedigree),file(params.bed))
		ZIPIT(params,burden_ch.zip.collect())
		}

workflow BURDEN_PAIRS {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()

		digest_ch = DIGEST_PAIRS(meta, reference, bed)
		version_ch = version_ch.mix(digest_ch.version)
	
		ch1_ch = BURDEN_SAMPLE_WGSELECT_PART_01(meta,reference,vcf, pedigree, digest_ch.merged_bed)
		version_ch = version_ch.mix(ch1_ch.version)


		each_gene = digest_ch.genes.splitText().map{it.trim()}
		

		header_ch = RVTESTS_REHEADER_01(meta, reference)
		version_ch = version_ch.mix(header_ch.version)

		assoc_ch = RVTEST_BY_PAIR(meta,
				reference,
				ch1_ch.contig_vcfs,
				ch1_ch.rvtest_pedigree,
				header_ch.output, 
				digest_ch.bed,
				each_gene.combine(each_gene).filter{T->T[0].compareTo(T[1])<=0}
				)
		version_ch = version_ch.mix(assoc_ch.version)
		
		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference, file(vcf).name ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden 1st intron", "Burden 1st intron ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip
	}

process DIGEST_PAIRS {
executor "local"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
output:
	path("merged.bed"),emit:merged_bed
	path("digest.bed"),emit:bed
	path("genes.txt"),emit:genes
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bedtools")}
set -o pipefail


mkdir -p TMP

${bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${bed}" | cut -f1-4 |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/digest.bed

# test no empty name
test `awk -F '\t' '(\$4=="")' TMP/digest.bed |wc -l` -eq 0

# extract uniq names
cut -f4 TMP/digest.bed | uniq | sort -T . | uniq > genes.txt

# merge all regions for wgselect
cut -f1,2,3 TMP/digest.bed |\
	LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > merged.bed


mv TMP/digest.bed ./


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">digest input bed, create gene list</entry>
	<entry key="versions">${getVersionCmd("bedtools jvarkit/bedrenamechr")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}

process RVTEST_BY_PAIR {
tag "${gene1} vs ${gene2}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(vcfs)
	path(pedigree)
	path(reheader)
	path(digest_bed)
	tuple val(gene1),val(gene2)
output:
	path("assoc.list"),emit:output
	path("version.xml"),emit:version
script:
	def  rvtest_params = "--burden 'cmc,exactCMC,zeggini' --kernel 'skato'"

	// --burden cmc,zeggini,mb,fp,exactCMC,cmcWald,rarecover,cmat --vt price,analytic --kernel 'skat[nPerm=1000],kbac,skato'

"""
hostname 1>&2
${moduleLoad("rvtests jvarkit bedtools bcftools")}
set -o pipefail
set -x
mkdir -p TMP ASSOC


awk -F '\t' '(\$4=="${gene1}" || \$4=="${gene2}")'  '${digest_bed}' |\
	cut -f1,2,3 |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/jeter.bed

test -s TMP/jeter.bed

bcftools concat -a --regions-file TMP/jeter.bed --file-list "${vcfs}" -O z -o TMP/jeter.vcf.gz
bcftools index --tbi TMP/jeter.vcf.gz

awk 'BEGIN {printf("${gene1}_${gene2}\t");} {printf("%s%s:%d-%s",(NR==1?"":","),\$1,int(\$2)+1,\$3);} END{printf("\\n");}' TMP/jeter.bed > TMP/jeter.setfile
test -s TMP/jeter.setfile


rvtest  --noweb \
        --inVcf TMP/jeter.vcf.gz \
	--setFile TMP/jeter.setfile \
	--pheno "${pedigree}" \
	--out "ASSOC/part." \
	${rvtest_params} 1>&2 2> TMP/last.rvtest.log

find \${PWD}/ASSOC -type f -name "part.*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for bed</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("bedtools rvtest bcftools")}</entry>
	<entry key="genes">${gene1}/${gene2}</entry>
</properties>
EOF
"""
}


process ZIPIT {
tag "N=${files.size()}"
publishDir "${meta.publishDir}" , mode: 'copy', overwrite: true
input:
	val(meta)
        val(files)
output:
        path("${meta.prefix?:""}archive.zip")
when:
        !meta.getOrDefault("publishDir","").trim().isEmpty()
script:
        prefix = meta.getOrDefault("prefix","")
"""

mkdir "${prefix}archive"

cat << EOF | while read F ; do ln -s "\${F}" "./${prefix}archive/" ; done
${files.join("\n")}
EOF

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}


runOnComplete(workflow);
