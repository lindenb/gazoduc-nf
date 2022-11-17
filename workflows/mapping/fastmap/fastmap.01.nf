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
params.fastqs = ""
params.help = false
params.publishDir = ""
params.prefix = ""
params.fastq_nchunks = 100

include {MAP_BWA_01} from '../../../subworkflows/mapping/map.bwa.01.nf'
include {runOnComplete;moduleLoad;getVersionCmd} from '../../../modules/utils/functions.nf'
include {COMPILE_FASTQ_SPLITTER} from '../../../modules/fastq/fastq.splitter.01.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'


def helpMessage() {
  log.info"""
## About

map fastqs on a reference genome

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --fastq_nchunks (int) how to split the fastq
  * --fastqs (file) one file containing the paths to the BAM/CRAM. Header: 'sample(tab)R1(tab)R2' [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume map.fastqs.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/referenceX.fasta \\
	--fastqs /path/to/input.tsv
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
	ch1_ch  = FASTMAP(params, params.reference,Channel.fromPath(params.fastqs))
	}


runOnComplete(workflow);


workflow FASTMAP {
	take:
		meta
		reference
		fastqs
	main:
		version_ch = Channel.empty()

		fastq_ch = fastqs.splitCsv(header:true,sep:'\t').
			flatMap{T->{
				int n = (params.fastq_nchunks as int);
				def array = [];
				for(int i=0;i< n;i++) {
					def h = T.plus("chunck":i,"n_chunks":n);
					array.add(h);	
					}
				return array;
				}}
		splitter_ch = COMPILE_FASTQ_SPLITTER(meta)
		version_ch = version_ch.mix(splitter_ch.version)

		bwa_ch = BWA(meta, reference, splitter_ch.executable, fastq_ch)
		version_ch = version_ch.mix(bwa_ch.version)

		merge_ch = SAMTOOLS_MERGE(meta, reference, bwa_ch.bam.groupTuple())
		version_ch = version_ch.mix(merge_ch.version)

		col_ch = COLLECT_TO_FILE_01([:],merge_ch.bam.map{T->T[1]}.collect())
		version_ch = version_ch.mix(col_ch.version)


		version_ch = MERGE_VERSION(meta, "fastmap", "fastmap",version_ch.collect())
		html = VERSION_TO_HTML(meta,version_ch)
	emit:
		version = version_ch
		html=html.html
		bams = col_ch.output
	}


process BWA {
tag "${row.sample} ${row.chunck}/${row.n_chunks} ${file(row.R1).name} ${row.R2?file(row.R2).name:""}"
afterScript "rm -rf TMP"
memory "5g"
cpus 5
input:
	val(meta)
	val(reference)
	val(splitter)
	val(row)
output:
	tuple val("${row.sample}"),path("${row.sample}.${row.chunck}.sorted.bam"),emit:bam
	path("version.xml"),emit:version
script:
	def sample = row.sample
	def CN = "Nantes"
	def PL = "ILLUMINA"
"""
hostname 1>&2
${moduleLoad("samtools/1.15.1 bwa")}
set -o pipefail
mkdir TMP

echo "\${LD_LIBRARY_PATH}" 1>&2 

${splitter} -n "${row.n_chunks}" -m "${row.chunck}" ${row.R1} ${row.R2?:""} |\
bwa mem -p -t ${task.cpus} -R '@RG\\tID:${sample}.${row.chunck}\\tSM:${sample}\\tLB:${sample}\\tCN:${CN}\\tPL:${PL}' "${reference}" - |\
	samtools view -O BAM -o TMP/jeter.bam


samtools sort --threads ${task.cpus} -o TMP/jeter2.bam -O BAM -T TMP/tmp TMP/jeter.bam
mv TMP/jeter2.bam TMP/jeter.bam


mv "TMP/jeter.bam" "${sample}.${row.chunck}.sorted.bam"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">map fastq to reference with bwa</entry>
	<entry key="sample">${sample}</entry>
	<entry key="reference">${reference}</entry>
	<entry key="bwa.version">${getVersionCmd("samtools bwa")}</entry>
</properties>
EOF
"""
}

process SAMTOOLS_MERGE {
tag "${sample} N=${L.size()}"
afterScript "rm -rf TMP"
memory "5g"
cpus 5
input:
	val(meta)
	val(reference)
	tuple val(sample),val(L)
output:
	tuple val(sample),path("${meta.prefix?:""}${sample}.merged.bam"),emit:bam
	path("version.xml"),emit:version
script:
"""
hostname 2>&1
${moduleLoad("samtools")}
set -o pipefail
mkdir TMP

samtools merge --reference "${reference}" --threads ${task.cpus} -o TMP/jeter.bam ${L.join(" ")}
samtools index -@ ${task.cpus} TMP/jeter.bam


mv TMP/jeter.bam "${meta.prefix?:""}${sample}.merged.bam"
mv TMP/jeter.bam.bai "${meta.prefix?:""}${sample}.merged.bam.bai" 

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">convert BAM to fastq</entry>
	<entry key="sample">${sample}</entry>
	<entry key="count(bam)">${L.size()}</entry>
        <entry key="versions">${getVersionCmd("samtools")}</entry>
</properties>
EOF
"""
}

