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

include {isHg19;runOnComplete;moduleLoad;getKeyValue;hasFeature;getVersionCmd;parseBoolean} from '../../../modules/utils/functions.nf'
include {BURDEN_SAMPLE_WGSELECT_PART_01}  from '../../../subworkflows/burden/burden.samples.wgselect.part.nf'
include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {RVTESTS_REHEADER_01} from '../../../modules/rvtests/rvtests.reheader.01.nf'
include {RVTESTS_POST_PROCESS} from '../../../subworkflows/rvtests/rvtests.post.process.01.nf'
include {CONCAT_FILES_01} from '../../../modules/utils/concat.files.nf'
include {SCATTER_TO_BED} from '../../../subworkflows/picard/picard.scatter2bed.nf'


params.reference=""
params.pedigree=""
params.vcf=""
params.disableFeatures="";
params.help=false
/** use only windows from this bed */
params.bed= "NO_FILE"
params.bed_cluster_method = " --size 1mb "
params.soacn = "" /* empty, take all */
/** params for bedtools makewindows */
params.makewindows_args = "-w 10000 -s 5000"

if(params.help) {
  log.info"""
## About

Burden for sliding windows

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
        --pedigree /path/to/input.ped \
```

## Workflow

![workflow](./workflow.svg)

"""
exit 0
}

workflow {
		burden_ch = BURDEN_SLIDING_WINDOWS(params, params.reference, params.vcf, file(params.pedigree),file(params.bed))
		ZIPIT(params,burden_ch.zip.collect())
		}

workflow BURDEN_SLIDING_WINDOWS {
	take:
		meta
		reference
		vcf
		pedigree
		bed
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()
		
		vcfbed_ch = VCF_TO_BED(meta,vcf)
		version_ch = version_ch.mix(vcfbed_ch.version)

		if(bed.name.equals("NO_FILE")) {	
			scatter_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"1000"],reference)
			version_ch = version_ch.mix(scatter_ch.version)
			utr_ch = EXTRACT_WINDOWS(meta, reference, vcfbed_ch.chromosomes, scatter_ch.bed)
			}
		else
			{
			utr_ch = EXTRACT_WINDOWS(meta, reference, vcfbed_ch.chromosomes, bed)
			}

		version_ch = version_ch.mix(utr_ch.version)
	
		ch1_ch = BURDEN_SAMPLE_WGSELECT_PART_01(meta,reference,vcf, pedigree, utr_ch.merged_bed)
		version_ch = version_ch.mix(ch1_ch.version)

		each_setfilelist_ch = utr_ch.output.splitText().
				map{S->file(S.trim())}
		
		header_ch = RVTESTS_REHEADER_01(meta, reference)
		version_ch = version_ch.mix(header_ch.version)

		assoc_ch = RVTEST_SLIDING(meta, reference, ch1_ch.contig_vcfs, ch1_ch.rvtest_pedigree, header_ch.output, each_setfilelist_ch)
		version_ch = version_ch.mix(assoc_ch.version)
		
		concat_ch = CONCAT_FILES_01(meta,assoc_ch.output.collect())
		version_ch = version_ch.mix(concat_ch.version)

		digest_ch = RVTESTS_POST_PROCESS(meta, reference, file(vcf).name ,concat_ch.output)
                version_ch = version_ch.mix(digest_ch.version)
		to_zip = to_zip.mix(digest_ch.zip)

		version_ch = MERGE_VERSION(meta, "burden UTR", "Burden UTR ${vcf}", version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params,version_ch.version)
		to_zip = to_zip.mix(html.html)
		
	emit:
		version = version_ch
		zip = to_zip
	}


process EXTRACT_WINDOWS {
memory "10g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(chromosomes)
	path(bed)
output:
	path("merged.bed"),emit:merged_bed
	path("setfile.list"),emit:output
	path("version.xml"),emit:version
script:
	def makewindows_args = meta.makewindows_args?:"-w 100000 -s 50000"
"""
hostname 1>&2
${moduleLoad("java bedtools htslib")}
set -o pipefail


mkdir -p TMP

# sort chromosomes
sort "${chromosomes}" > TMP/contigs.txt

# create merged bed
cut -f 1,2,3 '${bed}' |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	join -t '\t' -o '1.1,1.2,1.3' -1 1 -2 1 - TMP/contigs.txt |\
	bedtools merge > TMP/jeter.bed

# make windows
bedtools makewindows -b TMP/jeter.bed ${makewindows_args} |\
	awk -F '\t' '{printf("%s_%d_%s\t%s:%d-%s\\n",\$1,int(\$2)+1,\$3,\$1,int(\$2)+1,\$3);}' > TMP/jeter.setfile

# group the set files by pool of 'x' bases
mkdir BEDS
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" cluster --out BEDS  --size "1Mb" TMP/jeter.setfile

find \${PWD}/BEDS -type f -name "*.setfile" | sort -V -T TMP > setfile.list
test -s setfile.list

# make windows of 1Mb for wgselect
bedtools makewindows -w 1000000  -b TMP/jeter.bed > merged.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">make sliding windows</entry>
	<entry key="params">${makewindows_args}</entry>
	<entry key="versions">${getVersionCmd("bedtools tabix jvarkit/bedrenamechr jvarkit/setfiletools")}</entry>
	<entry key="bed">${bed}</entry>
</properties>
EOF
"""
}


process RVTEST_SLIDING {
tag "${setfile.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(vcfs)
	path(pedigree)
	path(reheader)
	path(setfile)
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


java -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" tobed "${setfile}" |\
	cut -f 1,2,3 | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.bed

#remove chr prefix
java -jar \${JVARKIT_DIST}/setfiletools.jar -R "${reference}" view --trim-chr '${setfile}' > TMP/jeter.setfile

bcftools concat -a --regions-file TMP/jeter.bed --file-list "${vcfs}" -O u |\
	bcftools annotate  --rename-chrs "${reheader}" -O z -o TMP/jeter.vcf.gz
bcftools index -t TMP/jeter.vcf.gz


rvtest  --noweb \
        --inVcf TMP/jeter.vcf.gz \
	--setFile TMP/jeter.setfile \
	--pheno "${pedigree}" \
	--out "ASSOC/part" \
	${rvtest_params} 1>&2 2> TMP/last.rvtest.log


find \${PWD}/ASSOC -type f -name "part*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for setfile</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("bedtools rvtest bcftools jvarkit/setfiletools")}</entry>
	<entry key="setfile">${setfile}</entry>
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
