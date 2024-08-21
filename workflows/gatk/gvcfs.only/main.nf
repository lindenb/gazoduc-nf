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

	beds_ch = Channel.fromPath(params.beds).
		splitText().
		map{file(it.trim())}

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
	hc_ch = HC_BAM_BED( genome_ch, dbsnp_ch, beds_ch.combine(bams_ch))
	if(params.xz_compression) {
		gather_ch = MERGE_VCFS( hc_ch.groupTuple().map{[it[0],it[1].plus(it[2])]})
		}
	else
		{
		gather_ch = XZ_VCFS( hc_ch.groupTuple().map{[it[0],it[1].plus(it[2])]})
		}
	}

process HC_BAM_BED {
tag "${bed.name} ${sample}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(dbsnp),path(dbsnp_tbi)
	tuple path(bed),val(sample),path(bam),path(bai)
output:
	tuple val(sample),path("*.g.vcf.gz"),path("*.g.vcf.gz.tbi"),emit:output
script:
	def prefix = bed.name.md5()+"."+sample;
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -Dsamjdk.compression_level=9" HaplotypeCaller \\
     -L "${bed}" \\
     -R "${fasta}" \\
     -I "${bam}" \\
     -ERC GVCF \\
     ${dbsnp.name.equals("NO_DBSNP")?"":"--dbsnp ${dbsnp}"} \\
     ${(params.mapq as Integer)<1?"":" --minimum-mapping-quality "+params.mapq} \\
     -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \\
     -O "TMP/jeter.g.vcf.gz"

mv TMP/jeter.g.vcf.gz "${prefix}.g.vcf.gz"
mv TMP/jeter.g.vcf.gz.tbi "${prefix}.g.vcf.gz.tbi"
"""
}

process MERGE_VCFS {
tag "${sample}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2

input:
	tuple val(sample),path("VCFS/*")
output:
	tuple path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi"),emit:output
script:
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

# cannot use GathVcfs because jvarkit cluster bed migh mix intervals
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -Dsamjdk.compression_level=9" MergeVcfs \\
	--COMPRESSION_LEVEL 9 \\
	--INPUT TMP/jeter.list \\
	--CREATE_INDEX true \\
	--OUTPUT TMP/jeter.g.vcf.gz

mv TMP/jeter.g.vcf.gz "${sample}.g.vcf.gz"
mv TMP/jeter.g.vcf.gz.tbi "${sample}.g.vcf.gz.tbi"
"""
}

process XZ_VCFS {
tag "${sample}"
label "process_quick"
afterScript "rm -rf TMP"
errorStrategy "retry"
maxRetries 2

input:
	tuple val(sample),path("VCFS/*")
output:
	tuple path("${sample}.g.vcf.xz"),emit:output
script:
"""
hostname 1>&2
module load gatk/0.0.0
mkdir -p TMP

find ./VCFS -name "*.vcf.gz" > TMP/jeter.list
test -s TMP/jeter.list

# cannot use GathVcfs because jvarkit cluster bed migh mix intervals
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" MergeVcfs \\
	--COMPRESSION_LEVEL 9 \\
	--INPUT TMP/jeter.list \\
	--CREATE_INDEX false \\
	--OUTPUT TMP/jeter.g.vcf

xz --force --best --extreme TMP/jeter.g.vcf
mv TMP/jeter.g.vcf.xz "${sample}.g.vcf.xz"
"""
}
