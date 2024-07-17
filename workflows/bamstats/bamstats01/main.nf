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

include {runOnComplete;moduleLoad} from '../../../modules/utils/functions.nf'


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


include {MOSDEPTH_DOWNLOAD_01} from '../../../modules/mosdepth/mosdepth.downoad.01.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'
include {SOMALIER_BAMS_01} from '../../../subworkflows/somalier/somalier.bams.01.nf'




workflow {
	version_ch = Channel.empty()
	toqc_ch = Channel.empty()
	tozip_ch = Channel.empty()
	bed_ch = Channel.fromPath(params.bed)
		


        samples0_ch = Channel.fromPath(params.samplesheet).
			splitCsv(header:true,sep:',').
			map{[it.sample,file(it.bam),file(it.bai)]}
	samples_ch = samples0_ch.combine(bed_ch)
	

        genome_ch = Channel.value([file(params.fasta),file(params.fasta+".fai"),file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict")])


	if(params.with_somalier) {


		som_ch = SOMALIER_BAMS_01(
			genome_ch,
			samples0_ch,
			file("NO_FILE") /* pedigree */,
			file("NO_FILE") /* custom sites */
			)
                //version_ch = version_ch.mix(som_ch.version)
		toqc_ch =  toqc_ch.mix(som_ch.qc)
		}





	if(params.with_samtools_idxstats) {
		fstats_ch = APPLY_SAMTOOLS_IDXSTAT(genome_ch,samples_ch)
		toqc_ch =  toqc_ch.mix(fstats_ch.output)
		}
	

	if(params.with_samtools_flagstat) {
		fstats_ch = APPLY_SAMTOOLS_FLAGSTAT(genome_ch,samples_ch)
		toqc_ch =  toqc_ch.mix(fstats_ch.output)
		}

	if(params.with_samtools_stats) {
		ststats_ch = APPLY_SAMTOOLS_STATS(genome_ch,samples_ch)
		toqc_ch =  toqc_ch.mix(ststats_ch.output)
		tozip_ch = tozip_ch.mix(ststats_ch.output.map{["samtools_stats",it]})
		}
 
	if(params.with_CollectWgsMetrics) {
		ch1_ch = APPLY_COLLECT_WGS_METRICS(genome_ch,samples_ch)
		version_ch = version_ch.mix(ch1_ch.version)
		toqc_ch =  toqc_ch.mix(ch1_ch.output)
		}

	if(params.with_mosdepth) {
                mosdepth_ch = MOSDEPTH_DOWNLOAD_01()
                version_ch = version_ch.mix(mosdepth_ch.version)
		
		
                ch2 = APPLY_MOSDEPTH(genome_ch,mosdepth_ch.executable,samples_ch)
                //version_ch = version_ch.mix(ch2.version)
		
		tozip_ch = tozip_ch.mix(ststats_ch.output.map{["mosdepth",it]})
		toqc_ch =  toqc_ch.mix(ch2.output)
		}
	
	ZIP_IT(tozip_ch.groupTuple())
         multiqc01_ch = MULTIQC_1(toqc_ch.collect())

         MULTIQC_2(file(params.sample2population), multiqc01_ch.json )
	
}



process APPLY_SAMTOOLS_IDXSTAT {
tag "${sample}"
label "process_quick"
cpus 1
input:
      	tuple path(fasta),path(fai),path(dict)
	tuple val(sample),path(bam),path(idx),path(bed)
output:
       	path("${sample}.idxstat.txt"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools")}

samtools idxstats --threads ${task.cpus} "${bam}" > "${sample}.idxstat.txt"
"""
}



process APPLY_SAMTOOLS_FLAGSTAT {
label "process_quick"
tag "${sample}"
afterScript "rm -rf TMP"
cpus 1
input:
      	tuple path(fasta),path(fai),path(dict)
	tuple val(sample),path(bam),path(idx),path(bed)
output:
       	path("${sample}.flagstats.txt"),emit:output
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail
mkdir -p TMP

samtools view ${bed.name.equals("NO_FILE")?"":"-M --regions-file \"${bed}\""} -O BAM --uncompressed --reference "${fasta}" "${bam}" |\
samtools flagstats '-' > TMP/jeter.flags.txt

mv -v "TMP/jeter.flags.txt" "${sample}.flagstats.txt" 

touch version.xml
"""
}


process APPLY_SAMTOOLS_STATS {
tag "${sample}"
label "process_quick"
time "3h"
afterScript "rm -rf TMP"
cpus 1
input:
      	tuple path(fasta),path(fai),path(dict)
	tuple val(sample),path(bam),path(idx),path(bed)
output:
       	path("${sample}.stats.txt"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail
mkdir -p TMP

samtools view ${bed.name.equals("NO_FILE")?"":"-M --regions-file \"${bed}\""}  -O BAM --uncompressed --reference "${fasta}" "${bam}" |\
samtools stats --ref-seq "${fasta}" --reference "${fasta}" > TMP/jeter.stats.txt

mv -v "TMP/jeter.stats.txt" "${sample}.stats.txt" 	

"""
}



process APPLY_COLLECT_WGS_METRICS {
tag "${sample}"
label "process_quick"
time "3h"
afterScript "rm -rf TMP"
memory "3g"
input:
      	tuple path(fasta),path(fai),path(dict)
	tuple val(sample),path(bam),path(idx),path(bed)
output:
        tuple path("${sample}.collect_wgs_metrics.txt"),emit:output
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0")}
set -o pipefail
mkdir -p TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  CollectWgsMetrics \
	I=${bam} \
	O=TMP/jeter.txt \
	R=${fasta}

mv -v TMP/jeter.txt "${sample}.collect_wgs_metrics.txt"
 

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">gatk CollectWgsMetrics</entry>
</properties>
EOF
"""
}


process APPLY_MOSDEPTH {
	tag "${sample} ${bam.name}"
	label "process_quick"
	time "3h"
	afterScript "rm -rf TMP"
	/** cpus 1 must be specified in config */
	input:
	      	tuple path(fasta),path(fai),path(dict)
		path(executable)
		tuple val(sample),path(bam),path(idx),path(bed)
	output:
		path("${sample}.*"),emit:output
	script:
		def mapq = 30
	"""
	hostname 1>&2
	mkdir -p TMP

        ./${executable} ${bed.name.equals("NO_FILE")?"":" --by \"${bed}\""} \\
 		-t ${task.cpus} --fasta "${fasta}" --mapq ${mapq} \\
		--use-median --no-per-base \\
		'TMP/${sample}' "${bam}"

	mv -v TMP/${sample}* ./

	"""
	}

process ZIP_IT {
	input:
		tuple val(name),path("FILES/*")
	output:
		path("${name}.zip"),emit:output
	script:
	"""
	zip -9j "${name}.zip" FILES/*
	"""
	}

process MULTIQC_1 {
label "process_medium"
input:
        path("FILES/*")
output:
        path("multiqc.01.zip"),emit:output
	path("${params.prefix?:""}multiqc/${params.prefix?:""}multiqc_report_data/multiqc_data.json"),emit:json
script:
        def prefix = params.prefix?:""
"""
hostname 1>&2
module load multiqc
mkdir -p TMP
find ./FILES -type l  > TMP/jeter.list

export LC_ALL=en_US.utf8

        mkdir -p "${prefix}multiqc"
        multiqc  --filename  "${prefix}multiqc_report.html" --no-ansi \
                        --title "BAM QCs"  \
                        --comment "QC for multiple BAMS"  \
                        --force \
                        --outdir "${prefix}multiqc" \
                        --file-list TMP/jeter.list

        rm -f multiqc.01.zip
        zip -9 -r "multiqc.01.zip" "${prefix}multiqc"

"""
}


process MULTIQC_2 {
label "process_medium"
input:
        path(sample2pop)
        path(json)
output:
        path("multiqc.02.zip"),emit:output
when:
	!sample2pop.name.equals("NO_FILE")
script:
        def prefix = params.prefix
"""
hostname 1>&2
module load jvarkit multiqc
mkdir -p TMP/DATA TMP/OUT2
cp -v '${json}' TMP/DATA/


export LC_ALL=en_US.utf8


java -jar \${JVARKIT_DIST}/jvarkit.jar multiqcpostproc --sample2collection "${sample2pop}" -o TMP/OUT2 TMP/DATA

find TMP/OUT2 -type f -name "*.json" > TMP/jeter.list


        mkdir -p "${prefix}multiqc.per.pop"
        multiqc  --filename  "${prefix}multiqc_report.html" --no-ansi \
                        --title "QC per population"  \
                        --comment "QC per population"  \
                        --force \
                        --outdir "${prefix}multiqc.per.pop" \
                        --file-list TMP/jeter.list

        rm -f multiqc.zip
        zip -9 -r "multiqc.02.zip" "${prefix}multiqc.per.pop"

"""
}
