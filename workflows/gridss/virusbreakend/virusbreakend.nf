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
params.bams=""
params.n_bam_per_call=1
params.help = false
params.conda = ""
/** publish Directory */
params.publishDir = ""
/** files prefix */
params.prefix = ""

include {isHg19;isHg38} from  '../../../modules/utils/functions.nf'
include {GRIDSS_SETUP_REFERENCE} from  '../../../modules/gridss/gridss.setupreference.nf'
include {SAMTOOLS_SAMPLES_01} from '../../../subworkflows/samtools/samtools.samples.01.nf'
include {SEQ_CACHE_POPULATE_01} from '../../../modules/samtools/seq.cache.populate.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {MERGE_VERSION} from '../../../modules/version/version.merge.nf'
include {SIMPLE_ZIP_01} from '../../../modules/utils/zip.simple.01.nf'

def helpMessage() {
  log.info"""
## About


## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list 
```

## Workflow

![workflow](./workflow.svg)
  
## See also

* 

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}

workflow {
	publish_ch = Channel.empty()

	gridss_ch = GRIDSS_VIRUS_BND(params,file(params.reference),Channel.fromPath(params.bams))
	}

workflow GRIDSS_VIRUS_BND {
	take:
		meta
		reference
		bams
	main:
		version_ch = Channel.empty()
		setup_ch = GRIDSS_SETUP_REFERENCE(meta,reference)
		version_ch = version_ch.mix(setup_ch.version)
		
		virus_ch = DOWNLOAD_VIRUS_BND(meta)
		version_ch = version_ch.mix(virus_ch.version)

		all_samples_ch = SAMTOOLS_SAMPLES_01(meta.plus("with_header":true), reference, file("NO_FILE"), bams)
		version_ch= version_ch.mix(all_samples_ch.version)

		each_sample_bam_ch = all_samples_ch.output.splitCsv(header:true,sep:'\t')

		cache_ch = SEQ_CACHE_POPULATE_01(meta, setup_ch.preproc_reference)
		version_ch= version_ch.mix(cache_ch.version)


		to_zip = Channel.empty()

		call_ch = CALL_VIRUS_BND(meta, setup_ch.preproc_reference, cache_ch.cache, virus_ch.virusdb, each_sample_bam_ch)
		version_ch = version_ch.mix(call_ch.version)
		to_zip = to_zip.mix(call_ch.vcf)

		version_ch = MERGE_VERSION(meta, "gridssvirusbnd", "GRIDSS Virus BND",version_ch.collect())
		to_zip = to_zip.mix(version_ch)

		html = VERSION_TO_HTML(meta,version_ch.version)
		to_zip = to_zip.mix(html.html)

		zip_ch = SIMPLE_ZIP_01(meta,to_zip.collect())
		
	emit:
		version = version_ch
		zip = zip_ch.zip
	}

process DOWNLOAD_VIRUS_BND {
afterScript "rm -f jeter.tar.gz"
input:
	val(meta)
output:
	path("virusbreakenddb_20210401"),emit:virusdb
	path("version.xml"),emit:version
script:
	def url = "https://virusbreakend.s3.us-east-2.amazonaws.com/virusbreakenddb_20210401.tar.gz"
script:
"""
wget -O jeter.tar.gz "${url}"
tar xvfz jeter.tar.gz

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">download Virus breakend data</entry>
       	<entry key="url"><a>${url}</a></entry>
</properties>
EOF
"""
}

process CALL_VIRUS_BND {
tag "${row.sample}"
afterScript "rm -rf TMP"
conda "${meta.conda}/gridss"
cpus 2
memory "10g"
input:
	val(meta)
	val(reference)
	val(refcache)
	path(virusdb)
	val(row)
output:
	path("${meta.prefix?:""}${row.sample}.output.vcf.gz"),emit:vcf
	path("version.xml"),emit:version
script:
	def sample=row.sample
	def bam=row.bam
"""
hostname 1>&2
mkdir -p TMP
export REF_CACHE=${refcache}/%2s/%2s/%s

virusbreakend \
  -r "${reference}" \
  --threads ${task.cpus} \
  --workingdir TMP \
  --db "${virusdb.toRealPath()}" \
  -o output.vcf \
  "${bam}"


grep "^#CHROM" -m1 output.vcf | tr "\t" "." 1>&2

# no name in the VCF?!
awk -F '\t' '/^#CHROM/ {OFS="\t";if(\$10=="") \$10="${sample}"; print;next;} {print;}' output.vcf > output2.vcf
mv output2.vcf output.vcf

find TMP -type f 1>&2

bcftools sort -T TMP -O z -o "${meta.prefix?:""}${sample}.output.vcf.gz" output.vcf
bcftools index -t "${meta.prefix?:""}${sample}.output.vcf.gz"
rm output.vcf

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">CALL virusdb with gridss</entry>
       	<entry key="reference">${reference}</entry>
       	<entry key="bam">${bam}</entry>
       	<entry key="sample">${sample}</entry>
</properties>
EOF
"""
}
