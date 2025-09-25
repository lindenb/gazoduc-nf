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
//include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include {runOnComplete;testKeyExistsAndNotEmpty;assertKeyMatchRegex} from '../../modules/utils/functions.nf'
include {BAM_QC                                                    } from '../../subworkflows/bamqc'
include {SCATTER_TO_BED                                            } from '../../subworkflows/gatk/scatterintervals2bed'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
//validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)



workflow {
	versions_ch = Channel.empty()
	toqc_ch = Channel.empty()
	tozip_ch = Channel.empty()
		
    def hash_ref= [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
      ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ hash_ref, file(params.fasta)]
	def fai   = [ hash_ref, file(params.fai)]
	def dict  = [ hash_ref, file(params.dict)]
 
  if(params.bed==null) {
    SCATTER_TO_BED(hash_ref,fasta,fai,dict)
    versions_ch = versions_ch.mix(SCATTER_TO_BED.out.versions)
    bed = SCATTER_TO_BED.out.bed
  } else {
      bed = Channel.of([hash_ref, file(params.bed)])
  }


    bams_ch = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.filter{testKeyExistsAndNotEmpty(it,"bam")}
			.filter{testKeyExistsAndNotEmpty(it,"sample")}
			.map{
				if(it.containsKey("bai")) return it;
				if(it.bam.endsWith(".bam")) return it.plus("bai":it.bam+".bai");
				return it.plus("bai":it.bam+".crai");
				}
            .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
            .map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
			.map{
				def ret=[id:it.sample];
				if(it.containsKey("population")) ret=ret.plus("population":it.population);
				if(it.containsKey("sex")) ret=ret.plus("sex":it.sex);
				if(it.containsKey("father")) ret=ret.plus("father":it.father);
				if(it.containsKey("mother")) ret=ret.plus("mother":it.mother);
				if(it.containsKey("family")) ret=ret.plus("family":it.family);
				return [ret,file(it.bam),file(it.bai)];
				}
	
	BAM_QC(
		hash_ref,
		fasta,
		fai,
		dict,
		bed,
		bams_ch
		)
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
