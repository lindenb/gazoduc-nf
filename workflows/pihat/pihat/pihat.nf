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
include {PIHAT01} from '../../../subworkflows/pihat/pihat.01.nf'
nextflow.enable.dsl=2

params.reference = ""
params.vcf = ""
params.help = false
params.publishDir = ""
params.prefix = ""
params.samples="NO_FILE"

def helpMessage() {
  log.info"""
## About

Pihat on vcf

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf (file) the VCF file [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume pihat.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--vcf /path/to/file.vcf
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

	pihat_ch = PIHAT01(
		params,
		params.reference,
		file(params.vcf),
		Channel.fromPath(params.samples)
		)
	//PUBLISH(indexcov_ch.zip)
	}


process PUBLISH {
tag "${zip.name}"
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(zip)
output:
	path(zip)
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
echo "publishing ${zip}"
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

