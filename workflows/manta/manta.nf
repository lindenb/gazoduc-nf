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
/** one file containing the paths to the BAM/CRAM  */
params.bams = ""
/** display help */
params.help = false
/** publish Directory */
params.publishDir = ""
/** files prefix */
params.prefix = ""
params.survivor_merge_params="1000 2 1 1 0 30"

include {MANTA_SINGLE_SV01} from '../../subworkflows/manta/manta.single.01.nf' 
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

def helpMessage() {
  log.info"""
## About

Detects CNV/SV using manta.

## Author

Pierre Lindenbaum PhD. Institut du Thorax. 44000 Nantes. France.

## Options

  * --reference (fasta) indexed fasta reference [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM. [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""
  * --survivor_merge_params (string) or empty to disable survivor. [${params.survivor_merge_params}]

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume manta.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also

  * https://github.com/Illumina/manta

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	manta_ch = MANTA_SINGLE_SV01(params, params.reference, params.bams)


	html = VERSION_TO_HTML(params,manta_ch.version)	

	to_publish = Channel.empty()
	to_publish = to_publish.
			mix(manta_ch.version).
			mix(manta_ch.zip).
			mix(manta_ch.survivor_vcf).
			mix(manta_ch.survivor_vcf_index).
			mix(manta_ch.manta_files).
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
	path("*.zip")
	path("*.list")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
	prefix = params.prefix?:""
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

