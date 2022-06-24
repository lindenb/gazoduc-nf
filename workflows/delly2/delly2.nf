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

include {DELLY2_RESOURCES} from '../../subworkflows/delly2/delly2.resources.nf' 
include {DELLY2_SV} from '../../subworkflows/delly2/delly2.sv.nf' 
include {SAMTOOLS_CASES_CONTROLS_01} from '../../subworkflows/samtools/samtools.cases.controls.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'

def helpMessage() {
  log.info"""
## About

Detects CNV/SV using delly2

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --cases (file) one file containing the paths to the BAM/CRAM for cases [REQUIRED]
  * --controls (file) one file containing the paths to the BAM/CRAM for controls [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""
  * --delly2_version (string) default: "${params.delly2_version}"
  * --cnv (boolean) Shall we call CNV ? default: "${params.cnv}"
  * --bnd (boolean) Shall we output BND ? default: "${params.bnd}"

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
	delly_ch = DELLY2_SV(params, params.reference, params.cases, params.controls)


	html = VERSION_TO_HTML(params,delly_ch.version)	

	to_publish = Channel.empty()
	to_publish = to_publish.
			mix(delly_ch.sv_vcf).
			mix(delly_ch.sv_vcf_index).
			mix(delly_ch.cnv_vcf).
			mix(delly_ch.cnv_vcf_index).
			mix(delly_ch.version).
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

