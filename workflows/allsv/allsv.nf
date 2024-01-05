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

/** path to indexed fasta reference */
params.reference = ""
/** one file containing the paths to the BAM/CRAM for controls */
params.controls = ""
/** one file containing the paths to the BAM/CRAM for cases */
params.cases = ""
/** display help */
params.help = false
/** publish Directory */
params.publishDir = ""
/** files prefix */
params.prefix = ""
/** delly version *:
params.delly2_version = "v1.0.3"
/** keep INFO/TYPE=BND */
params.bnd=true
/** search CNVs */
params.cnv=true

include {DELLY2_SV} from '../../subworkflows/delly2/delly2.sv.nf' 
include {INDEXCOV} from '../../subworkflows/indexcov/indexcov.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'


def helpMessage() {
  log.info"""
## About

Detects SV using Delly2, manta, smoove, indexcov.

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --cases (file) one file containing the paths to the BAM/CRAM for cases [REQUIRED]
  * --controls (file) one file containing the paths to the BAM/CRAM for controls [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume ${workflow.scriptFile} \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--cases /path/to/bams.cases.list \\
	--controls /path/to/bams.controls.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also

  * https://github.com/dellytools/delly

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	version = Channel.emtpy()
	
	delly_ch = DELLY2_SV(params, params.reference, params.cases, params.controls)
	version_ch = version_ch.mix(delly_ch.version)

	smoove_ch = SMOOVE_SV_POPULATION_01(params, params.reference, params.cases, params.controls)
	version_ch = version_ch.mix(smoove_ch.version)

	

	html = VERSION_TO_HTML(params,smoove_ch.version)	


	concat_ch =Channel.empty()
	concat_ch = concat_ch.mix(params.cases)
	if(!isBlank(params.controls)) {
		concat_ch = concat_ch.mix(params.controls)
		}
	

	cat_files_ch = COLLECT_TO_FILE_01(params,concat_ch.collect())
	indexcov_ch = INDEXCOV(params,params.reference,cat_files_ch.output)
	version_ch = version_ch.mix(indexcov_ch.version)

	manta_ch = MANTA_SINGLE_SV01(params, params.reference, files_ch.output)
	version_ch = version_ch.mix(manta_ch.version)

	version_ch = MERGE_VERSION("allsv","several structural variation tools...",version_ch.collect())

	html = VERSION_TO_HTML(params,version_ch)	

	to_publish = Channel.empty()
	to_publish = to_publish.
			mix(delly_ch.sv_vcf).
			mix(delly_ch.sv_vcf_index).
			mix(delly_ch.cnv_vcf).
			mix(delly_ch.cnv_vcf_index).
			mix(version_ch).
			mix(html.html)

	PUBLISH(to_publish.collect())
	}

process PUBLISH {
tag "N=${files.size()}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(files)
output:
	path("*.bcf")
	path("*.bcf.csi")
	path("*.xml")
	path("*.html")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
	prefix = params.getOrDefault("prefix","")
"""
for F in ${files.join(" ")}
do
	ln -s "\${F}" ./
done
"""
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

