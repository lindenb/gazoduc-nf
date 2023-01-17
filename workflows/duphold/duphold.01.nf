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

/** path to indexed fasta reference */
params.reference = ""
params.references = ""
/* the multi sample vcf or a file with the .list suffix */
params.vcf = "NO_FILE"
params.bams = "NO_FILE"
/** optional path to snp/indel VCF/BCF with which to annotate SVs. BCF is highly recommended as it's much faster to parse. */
params.snps = "NO_FILE"
params.help = false
params.publishDir = ""
params.prefix = ""

include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;isHg19;isHg38;runOnComplete;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

def helpMessage() {
  log.info"""
## About

apply mosdepth to a set of bams.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --vcf (file) VCF file
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow1.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--vcf /path/to/vcf
```

## Workflow

![workflow](./workflow.svg)
  
## See also


"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	ch1 = DUPHOLD01(params,params.reference,file(params.bams),file(params.vcf),file(params.snps))
	}

runOnComplete(workflow)

workflow DUPHOLD01 {
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
	
		join_ch = JOIN_BAM_VCFS(meta,bams_ch.output,vcf)
		version_ch = version_ch.mix(join_ch.version)

		exec_ch = DOWNLOAD_DUPHOLD(meta)
		version_ch = version_ch.mix(exec_ch.version)

		per_contig =join_ch.output.splitCsv(header:true,sep:'\t').
			flatMap{T->T.contigs.split("[,]").collect{V->T.plus(contig:V)} }

		to_zip = Channel.empty()
		call_ch= APPLY_DUPHOLD(meta, exec_ch.executable, snps,per_contig )
		version_ch = version_ch.mix(call_ch.version)
		
		concat_ch = CONCAT(meta, call_ch.vcf.groupTuple())
		version_ch = version_ch.mix(concat_ch.version)
		to_zip = to_zip.mix(concat_ch.vcf)
		
		version_ch = MERGE_VERSION(meta, "Duphold", "Duphold",version_ch.collect())		
		to_zip = to_zip.mix(version_ch)

		html = VERSION_TO_HTML(meta,version_ch.version)
		to_zip = to_zip.mix(html.html)

		zip_ch = SIMPLE_ZIP_01(meta,to_zip.collect())
		
	emit:
		zip = zip_ch.zip
		version= version_ch
	}



process JOIN_BAM_VCFS {
tag "${vcf} ${bams}"
afterScript "rm -rf TMP"
input:
	val(meta)
	path(bams)
	path(vcf)
output:
	path("join.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir -p TMP OUT

# path to all the vcfs

if ${vcf.name.endsWith(".list")} ; then
        cp "${vcf}" TMP/jeter.list
else
	echo '${vcf.toRealPath()}' > TMP/jeter.list
fi

# sort BAMs by name
sort -t '\t' -k1,1 "${bams}" > TMP/jeter.a

echo -e "vcf\tid\tcontigs\treference\tbams" > join.tsv

test -s TMP/jeter.list
i=1
cat TMP/jeter.list | while read V
do
	bcftools query -l "\${V}" | sort > TMP/jeter.samples


	echo -e -n "\${V}\t" >> join.tsv


	# find a name for this vcf 
	if test `wc -l <  TMP/jeter.samples` -eq 1 ; then
		cat TMP/jeter.samples | tr "\\n" "\t" >> join.tsv
	else
		echo -e -n "\${i}\t" >> join.tsv
	fi
	
	# extract contigs
	bcftools index -s "\${V}" | cut -f1 | paste -s -d ',' | tr "\\n" "\t" >> join.tsv

	# print reference
	join -t '\t' -1 1 -2 1 -o '1.4' TMP/jeter.a TMP/jeter.samples | sort | uniq | tr "\\n" "\t" >> join.tsv
	# make list of bams
	join -t '\t' -1 1 -2 1 -o '1.3' TMP/jeter.a TMP/jeter.samples | sort | uniq > "OUT/bams.\${i}.list"
	test -s OUT/bams.\${i}.list
	echo "\${PWD}/OUT/bams.\${i}.list" >> join.tsv

	i=\$((i+1))
done


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join vcf and bams</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="bams">${bams}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
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
tag "contig:${row.contig} ID:${row.id} ${file(row.vcf).name} ${file(row.bams).name}"
afterScript "rm -rf TMP"
memory "10g"
cpus 4
input:
	val(meta)
	val(duphold)
	path(snps)//TODO use it
	val(row)
output:
	tuple val("${row.id}"),path("${meta.prefix?:""}${row.id}.${row.contig}.duphold.bcf"),emit:vcf
	path("version.xml"),emit:version
script:

"""
hostname 1>&2
${moduleLoad("bcftools")}
mkdir TMP

bcftools view --threads "${task.cpus}" -O b -o TMP/jeter.bcf "${row.vcf}" "${row.contig}"

cat "${row.bams}" | while read B
do

echo "#Processing \${B}" 1>&2

${duphold} --vcf "TMP/jeter.bcf" \
	--threads "${task.cpus}"  \
	--bam "\${B}" \
	--fasta "${row.reference}" \
	--output TMP/jeter2.bcf 1>&2

mv -v TMP/jeter2.bcf  TMP/jeter.bcf

done

bcftools view --threads "${task.cpus}" -O b -o "${meta.prefix?:""}${row.id}.${row.contig}.duphold.bcf" TMP/jeter.bcf
bcftools index --threads "${task.cpus}" "${meta.prefix?:""}${row.id}.${row.contig}.duphold.bcf"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">apply duphold</entry>
        <entry key="reference">${row.reference}</entry>
        <entry key="contig">${row.contig}</entry>
        <entry key="vcf">${row.vcf}</entry>
        <entry key="id">${row.id}</entry>
        <entry key="bams">${row.bams}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}

process CONCAT {
tag "${id} N=${L.size()}"
input:
	val(meta)
	tuple val(id),val(L)
output:
	path("${meta.prefix?:""}${id}.duphold.bcf"),emit:vcf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}


bcftools concat --allow-overlaps -O b -o "${meta.prefix?:""}${id}.duphold.bcf" ${L.join(" ")}
bcftools index "${meta.prefix?:""}${id}.duphold.bcf"

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">apply duphold</entry>
        <entry key="id">${id}</entry>
        <entry key="contigs.count">${L.size()}</entry>
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
