/*

Copyright (c) 2024 Pierre Lindenbaum

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

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)




workflow {

	bed_ch = Channel.fromPath(params.bed).
		splitCsv(sep:'\t',header:false).
		map{[it[0],(it[1] as int)+1,it[2]]}

	genome_ch = [file(params.fasta) , file(params.fasta+".fai"), file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict" ) ]	


	bams_ch = Channel.fromPath(params.samplesheet).
		splitCsv(header:true,sep:'\t').
		map{
			[
			it.sample,
			file(it.bam),
			file(it.bam_index?:it.bam+(it.bam.endsWith(".bam")?".bai":".crai"))
			]
		}.view()


	if(params.dbsnp.equals("NO_FILE"))
		{
		dbsnp_ch = [file("NO_DBSNP"), file("NO_DBSNP_TBI")]
		}
	else
		{
		dbsnp_ch = [file(params.dbsnp), file(params.dbsnp+".tbi")]
		}
	hc_ch = HC_BAM_BED( genome_ch, dbsnp_ch, bed_ch.combine(bams_ch))

	if(params.compression.equals("xz")) {
		gather_ch = XZ_VCFS( hc_ch.groupTuple().map{[it[0],it[1].plus(it[2])]})
		}
	else if(params.compression.equals("bz2") || params.compression.equals("bzip2")) {
		gather_ch = BZ2_VCFS( hc_ch.groupTuple().map{[it[0],it[1].plus(it[2])]})
		}
	else
		{
		gather_ch = MERGE_VCFS( hc_ch.groupTuple().map{[it[0],it[1].plus(it[2])]})
		}
	}

process HC_BAM_BED {
tag "${chrom}:${start}-${end} ${sample}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(dbsnp),path(dbsnp_tbi)
	tuple val(chrom),val(start),val(end),val(sample),path(bam),path(bai)
output:
	tuple val(sample),path("*.g.vcf.gz"),path("*.g.vcf.gz.tbi"),emit:output
script:
	def interval = "${chrom}:${start}-${end}"
	def prefix = interval.md5()+"."+sample;
	def arg1 = task.ext.arg1?:""
	def arg2 = task.ext.arg2?:""
	def compression_level = task.ext.compression_level?:5;
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

   gatk --java-options "${arg1} -XX:-UsePerfData -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -Dsamjdk.compression_level=${compression_level}" HaplotypeCaller \\
     -L "${interval}" \\
     -R "${fasta}" \\
     -I "${bam}" \\
     -ERC GVCF \\
     ${dbsnp.name.equals("NO_DBSNP")?"":"--dbsnp ${dbsnp}"} \\
     ${(params.mapq as Integer)<1?"":" --minimum-mapping-quality "+params.mapq} \\
     ${arg2} \\
     -O "TMP/jeter.g.vcf.gz"

mv TMP/jeter.g.vcf.gz "${prefix}.g.vcf.gz"
mv TMP/jeter.g.vcf.gz.tbi "${prefix}.g.vcf.gz.tbi"
"""
}

process GATHER_VCFS {
tag "${sample}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2

input:
	tuple val(sample),path("VCFS/*")
output:
	path("${sample}.g.vcf.gz"),emit:output
script:
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

gatk --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP  -Dsamjdk.compression_level=9" GatherVcfs \\
	--REORDER_INPUT_BY_FIRST_VARIANT \\
	--COMPRESSION_LEVEL 9 \\
	--INPUT TMP/jeter.list \\
	--OUTPUT TMP/jeter.g.vcf.gz

mv TMP/jeter.g.vcf.gz "${sample}.g.vcf.gz"
"""
}

process XZ_VCFS {
tag "${sample}"
label "process_medium"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2

input:
	tuple val(sample),path("VCFS/*")
output:
	path("${sample}.g.vcf.xz"),emit:output
script:
	def args1 = task.ext.args1?:"--force --best --extreme"
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

gatk --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP  -Dsamjdk.compression_level=9" GatherVcfs \\
	--REORDER_INPUT_BY_FIRST_VARIANT \\
	--INPUT TMP/jeter.list \\
	--OUTPUT TMP/jeter.g.vcf

xz -v --threads=${task.cpus} ${args1} TMP/jeter.g.vcf > xv.log
mv TMP/jeter.g.vcf.xz "${sample}.g.vcf.xz"
"""
}


process BZ2_VCFS {
tag "${sample}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2

input:
	tuple val(sample),path("VCFS/*")
output:
	path("${sample}.g.vcf.bz2"),emit:output
script:
	def args1 = task.ext.args1?:" --best "
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

gatk --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP  -Dsamjdk.compression_level=9" GatherVcfs \\
	--REORDER_INPUT_BY_FIRST_VARIANT \\
	--INPUT TMP/jeter.list \\
	--OUTPUT TMP/jeter.g.vcf

ls -lah TMP/jeter.g.vcf 1>&2

bzip2 -v  ${args1} TMP/jeter.g.vcf > bz2.log
mv TMP/jeter.g.vcf.bz2 "${sample}.g.vcf.bz2"
"""
}
