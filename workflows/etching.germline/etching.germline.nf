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
nextflow.enable.dsl=2

/** path to indexed fasta reference */
params.reference = ""
params.fastqs = ""
params.help = false
params.publishDir = ""
params.bed = "NO_FILE"
params.prefix = ""

include {moduleLoad;runOnComplete; getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

def helpMessage() {
  log.info"""
## About

run Etching for a list of fastqs

## Author

Pierre Lindenbaum

## Options

  * --fastqs (file) tab delimited [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--fastqs input.tsv
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
	ch1 = ETCHING_01(params, params.reference, file(params.bed), Channel.fromPath(params.fastqs))
	html = VERSION_TO_HTML(params, ch1.version)	
	//PUBLISH(ch1.zip)
	}

process PUBLISH {
executor "local"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(zip)
output:
	path("*.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s ${zip} ./
"""
}

runOnComplete(workflow);

workflow ETCHING_01 {
	take:
		meta
		reference
		bed
		fastqs
	main:
		version_ch = Channel.empty()

		kmers_ch = DOWNLOAD_KMERS([])
		version_ch = version_ch.mix(kmers_ch.version)

		etching_ch = COMPILE_ETCHING([])
		version_ch = version_ch.mix(etching_ch.version)


		fq_ch = fastqs.splitCsv(header:true,sep:',')

		sv_ch = APPLY_ETCHING(meta, reference, etching_ch.executable, kmers_ch.output, bed, fq_ch)
		version_ch = version_ch.mix(sv_ch.version)

		zip_ch = SIMPLE_ZIP_01([:] ,sv_ch.output.flatMap{T->[T[1],T[2]]}.collect())
		version_ch = version_ch.mix(zip_ch.version)
	
		version_ch = MERGE_VERSION(meta, "etching", "Etching", version_ch.collect())
	emit:
		version = version_ch
		zip = zip_ch.zip
	}


process DOWNLOAD_KMERS {
afterScript "rm -f PGK2.tar.gz"
errorStrategy "finish"
input:
	val(meta)
output:
	path("PGK2/PGK2"),emit:output
	path("version.xml"),emit:version
script:
	def url="http://big.hanyang.ac.kr/ETCHING/PGK2.tar.gz"
"""
hostname 1>&2

wget -O "PGK2.tar.gz" "${url}"
tar xvfz PGK2.tar.gz
# base file
touch PGK2/PGK2

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">download kmers</entry>
	<entry key="url">${url}</entry>
	<entry key="version">${getVersionCmd("wget")}</entry>
</properties>
EOF
"""
}

process COMPILE_ETCHING {
input:
	val(meta)
output:
	path("ETCHING/bin/etching"),emit:executable
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bwa samtools")}

git clone -b patch-1 --depth=1 "https://github.com/lindenb/ETCHING.git"
(cd ETCHING && make)


#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">compile etching</entry>
	<entry key="version">`./ETCHING/bin/etching --version | paste -sd ' '`</entry>
</properties>
EOF
"""
}

process APPLY_ETCHING {
tag "${row.sample}"
afterScript "rm -rf TMP"
cpus 10
memory "15g"
input:
	val(meta)
	val(reference)
	path(executable)
	path(kmers)
	path(bed)
	val(row)
output:
	tuple val(row), path("${meta.prefix?:""}${row.sample}.TF.BND.vcf.gz"), path("${meta.prefix?:""}${row.sample}.TF.SV.vcf.gz"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools bwa samtools")}
mkdir -p TMP

${executable.toRealPath()} \
		-1 "${row.R1}" -2 "${row.R2}" \
		-f "${kmers.toRealPath()}" \
		${bed.name.equals("NO_FILE")?"":"--target-filter ${bed}"} \
		-g "${reference}" -t ${task.cpus} -o "TMP"


find \${PWD}/ 1>&2

cat << EOF > TMP/rename.txt
FIRST_MATE	FIRST_MATE_${row.sample}
SECOND_MATE	SECOND_MATE_${row.sample}
EOF

bcftools reheader  --samples TMP/rename.txt TMP.TF.BND.vcf |\
bcftools sort -T TMP -O z -o "${meta.prefix?:""}${row.sample}.TF.BND.vcf.gz"

bcftools reheader  --samples TMP/rename.txt TMP.TF.SV.vcf |\
bcftools sort -T TMP -O z -o "${meta.prefix?:""}${row.sample}.TF.SV.vcf.gz"

#######################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">apply etching</entry>
	<entry key="sample">${row.sample}</entry>
	<entry key="R1">${row.R1}</entry>
	<entry key="R2">${row.R2}</entry>
	<entry key="reference">${reference}</entry>
	<entry key="bed">${bed}</entry>
	<entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}

