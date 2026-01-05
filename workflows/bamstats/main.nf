/*

Copyright (c) 2026 Pierre Lindenbaum

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
include {assertKeyExistsAndNotEmpty                                } from '../../modules/utils/functions.nf'
include {BAM_QC                                                    } from '../../subworkflows/bamqc'
include {SCATTER_TO_BED                                            } from '../../subworkflows/gatk/scatterintervals2bed'
include {MULTIQC                                                   } from '../../subworkflows/multiqc'
include {PREPARE_ONE_REFERENCE                                     } from '../../subworkflows/samtools/prepare.one.ref'
include {PREPARE_USER_BED                                          } from '../../subworkflows/bedtools/prepare.user.bed'

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
	versions = Channel.empty()
	multiqc = Channel.empty()

   def project_root = file("${launchDir}")
   if(params.fasta==null) {
	throw new IllegalArgumentException("--fasta undefined")
	}
   
   if(params.samplesheet==null) {
	throw new IllegalArgumentException("--samplesheet undefined")
	}
		
    def metadata = [
      id: "bamstats"
      ]
  
  
  PREPARE_ONE_REFERENCE(metadata, Channel.of(file(params.fasta)).map{f->[[id:f.baseName],f]})
  versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)

  if(params.bed==null) {
    bed = PREPARE_ONE_REFERENCE.out.scatter_bed
  } else {

	PREPARE_USER_BED(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			PREPARE_ONE_REFERENCE.out.scatter_bed,
			Channel.of([[id:file(params.bed).baseName],file(params.bed)])
			)
	versions = versions.mix(PREPARE_USER_BED.out.versions)
	bed = PREPARE_USER_BED.out.bed.first()
  	}


    bams_ch = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.map{assertKeyExistsAndNotEmpty(it,"bam")}
			.map{assertKeyExistsAndNotEmpty(it,"sample")}
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
				return [ret,project_root.resolve(file(it.bam)), project_root.resolve(file(it.bai))];
				}
	
	BAM_QC(
		metadata.plus(
			with_mosdepth: params.with_mosdepth,
			with_samtools_stats : params.with_samtools_stats,
			with_samtools_idxstats : params.with_samtools_idxstats,
			with_samtools_flagstats : params.with_samtools_flagstat,
			with_collect_metrics : params.with_CollectWgsMetrics,
			with_depth_outliers : params.with_depth_outliers
			),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		bed,
		bams_ch
		)
	versions = versions.mix(BAM_QC.out.versions)
	multiqc  = multiqc.mix(BAM_QC.out.multiqc)

	
	MULTIQC(
		metadata,
		Channel.of([[id:"no_sample2col"],[]]),
		versions,
		Channel.of([[id:"no_mqc_config"],[]]),
		multiqc
		)

}

runOnComplete(workflow)



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
