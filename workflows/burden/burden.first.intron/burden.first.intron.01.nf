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
params.disableFeatures="";
params.help=false
/** use only introns overlapping this bed */
params.bed= "NO_FILE"
/** in  intron, only keep/annotate in this interval(s) */
params.regulatory_bed = "NO_FILE2" //yes FILE2 to avoid conflicts for symbolic links
params.bed_cluster_method = " --size 1mb "
params.soacn = "SO:0001627,SO:0001568"

if(params.help) {
  log.info"""
## About

Burden for 1st intron.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
  * --bed <file> optional bed file to limit the analysis to the genes overlapping a  bed file.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume burden.first.intron.01.nf \\
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
		burden_ch = BURDEN_1ST_INTRON(params, params.reference, params.vcf, file(params.pedigree),file(params.bed),file(params.regulatory_bed))
		ZIPIT(params,burden_ch.zip.collect())
		}

workflow BURDEN_1ST_INTRON {
	take:
		meta
		reference
		vcf
		pedigree
		bed
		regulatory_bed
	main:

		version_ch = Channel.empty()
		to_zip = Channel.empty()

		introns_ch = EXTRACT_1ST_INTRON(meta, reference, bed, regulatory_bed)
		version_ch = version_ch.mix(introns_ch.version)
	
		ch1_ch = BURDEN_SAMPLE_WGSELECT_PART_01(meta,reference,vcf, pedigree, introns_ch.merged_bed)
		version_ch = version_ch.mix(ch1_ch.version)


		each_intron_bed_ch = introns_ch.output.splitText().map{it.trim()}.map{S->file(S)}
		
		header_ch = RVTESTS_REHEADER_01(meta, reference)
		version_ch = version_ch.mix(header_ch.version)

		assoc_ch = RVTEST_FIRST_INTRON(meta, reference, ch1_ch.contig_vcfs, ch1_ch.rvtest_pedigree, header_ch.output, each_intron_bed_ch)
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

process EXTRACT_1ST_INTRON {
memory "3g"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(bed)
	path(regulatory_bed)
output:
	path("merged.bed"),emit:merged_bed
	path("first.introns.list"),emit:output
	path("version.xml"),emit:version
script:
	def url = isHg19(reference)?"http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ensGene.txt.gz":""
"""
hostname 1>&2
${moduleLoad("jvarkit bedtools")}
set -o pipefail

test ! -z "${url}"

mkdir -p TMP

wget -O - "${url}" | gunzip -c |\
	awk -F '\t' '(int(\$9)>1) {split(\$10,S,/[,]/);split(\$11,E,/[,]/);N=int(\$9);if(\$4=="+") printf("%s\t%d\t%d\t%s\\n",\$3,E[1],int(S[2]),\$2);if(\$4=="-") printf("%s\t%d\t%d\t%s\\n",\$3,E[N-1],int(S[N]),\$2);}' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  |\
	sort -T TMP -t '\t' -k1,1 -k2,2n -k3,3n --unique > TMP/jeter2.bed


if ${!bed.name.equals("NO_FILE")} ; then

	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP  "${bed}" |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n | cut -f1,2,3 > TMP/jeter3.bed
	test -s TMP/jeter3.bed

	bedtools intersect -u -a TMP/jeter2.bed -b TMP/jeter3.bed |\
	LC_ALL=C sort -S ${task.memory.kilo} -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter4.bed
	mv TMP/jeter4.bed TMP/jeter2.bed

fi

# group the BED of each introns by pool of 'x' bases
mkdir BEDS
java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/bedcluster.jar -R "${reference}" --out BEDS --names --size "1Mb" TMP/jeter2.bed
find \${PWD}/BEDS -type f -name "*.bed" > first.introns.list
test -s first.introns.list

# create merged bed

cut -f 1,2,3 TMP/jeter2.bed |\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\
	bedtools merge > TMP/jeter.bed

if ${!regulatory_bed.name.equals("NO_FILE2") /* yes 'NO_FILE2' to avoid conflicts with syminks */ } ; then

	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP "${regulatory_bed}" |\
		cut -f1,2,3 |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge  > TMP/jeter2.bed

	bedtools intersect -a TMP/jeter.bed -b TMP/jeter2.bed |\
		LC_ALL=C sort -T TMP -t '\t' -k1,1 -k2,2n > TMP/jeter3.bed

	mv TMP/jeter3.bed TMP/jeter.bed	

fi

mv TMP/jeter.bed merged.bed

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">1st intron from ucsc</entry>
	<entry key="url"><url>${url}</url></entry>
	<entry key="versions">${getVersionCmd("wget bedtools jvarkit/bedrenamechr jvarkit/bedcluster")}</entry>
	<entry key="bed">${bed}</entry>
	<entry key="regulatory_bed">${regulatory_bed}</entry>
</properties>
EOF
"""
}

process RVTEST_FIRST_INTRON {
tag "${bed.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	path(vcfs)
	path(pedigree)
	path(reheader)
	path(bed)
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

cut -f 1,2,3 "${bed}" | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter.bed
bcftools concat -a --regions-file TMP/jeter.bed --file-list "${vcfs}" -O b -o TMP/jeter.bcf
bcftools index TMP/jeter.bcf

i=1
awk -F '\t' '{printf("%s:%d-%s\t%s\\n",\$1,int(\$2)+1,\$3,\$4);}' "${bed}" | while read RGN ENST
do
	bcftools view TMP/jeter.bcf "\${RGN}" |\
		java -jar \${JVARKIT_DIST}/vcfburdenfiltergenes.jar -a "\${ENST}" |\
		bcftools annotate -O z --rename-chrs "${reheader}" -o TMP/jeter.vcf.gz -

	if test `bcftools query -f '.' TMP/jeter.vcf.gz | wc -c` -gt 0 ; then

		bcftools index --tbi -f TMP/jeter.vcf.gz
	

		echo -n "## \${ENST}: " && bcftools query -f '.' TMP/jeter.vcf.gz | wc -c

		# build setFile
		echo "\${ENST}\t\${RGN}" | sed 's/\tchr/\t/' > TMP/variants.setfile

		rvtest  --noweb \
        		--inVcf TMP/jeter.vcf.gz \
			--setFile TMP/variants.setfile \
	        	--pheno "${pedigree}" \
		        --out "ASSOC/part.\${i}" \
			${rvtest_params} 1>&2 2> TMP/last.rvtest.log

		i=\$((i+1))

	fi
done

find \${PWD}/ASSOC -type f -name "part.*assoc" > assoc.list

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">invoke rvtest for bed</entry>
	<entry key="rvtest.path">\$(which rvtest)</entry>
	<entry key="versions">${getVersionCmd("bedtools rvtest bcftools jvarkit/vcfburdenfiltergenes")}</entry>
	<entry key="bed">${bed}</entry>
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
