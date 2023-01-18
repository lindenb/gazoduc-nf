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
params.smoove_image = ""
params.gff3 = ""
/* use duphold ? slow if too many samples */
params.with_duphold = false

include {SMOOVE_SV_POPULATION_01} from '../../subworkflows/smoove/smoove.population.01.nf' 
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


def helpMessage() {
  log.info"""
## About

Detects CNV/SV using smoove

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --cases (file) one file containing the paths to the BAM/CRAM for cases [REQUIRED]
  * --controls (file) one file containing the paths to the BAM/CRAM for controls
  * --publishDir (dir) Save output in this directory. default: "${params.publishDir}"
  * --prefix (string) files prefix. default: "${params.prefix}"
  * --smoove_image (string) path to precompiled singularuty image default: "${params.smoove_image}"
  * --gff3 (path/url) path to gff3 file for annotation. If blank, the default file is downloaded from the web.

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume ${file(workflow.scriptFile)} \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--cases /path/to/bams.cases.list \\
	--controls /path/to/bams.controls.list
```

## Workflow

![workflow](./workflow.svg)
  
## See also

  * https://github.com/brentp/smoove

"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	smoove_ch = SMOOVE_SV_POPULATION_01(
		params,
		params.reference,
		params.cases,
		params.controls
		)

	html = VERSION_TO_HTML(params,smoove_ch.version)	

	to_publish = Channel.empty()
	to_publish = to_publish.
			mix(smoove_ch.vcf).
			mix(smoove_ch.index).
			mix(smoove_ch.version).
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

runOnComplete(workflow)
