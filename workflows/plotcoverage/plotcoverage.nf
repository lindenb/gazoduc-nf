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

def gazoduc = gazoduc.Gazoduc.getInstance().putDefaults(params)


/** path to indexed fasta reference */
//params.reference = ""
params.references = "NO_FILE"
params.mapq = 0
params.bams = ""
params.bed = "BED"
//params.help = false
//params.publishDir = ""
//params.prefix = ""


include {PLOT_COVERAGE_01} from '../../subworkflows/plotdepth/plot.coverage.01.nf'
include {runOnComplete} from '../../modules/utils/functions.nf'


gazoduc.reference().required().put(params)

def helpMessage() {
  log.info"""
## About

apply mosdepth to a set of bams.

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --vcf (file) required SV indexed vcf file. default: ""
  * --excludeids (file) optional file containing IDs of SV to ignore.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume plotcoverage.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--bed x.bed
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

gazoduc.usage().print()

workflow {
	ch1 = PLOT_COVERAGE_01(params,params.reference, file(params.references),file(params.bams), file(params.bed))
	//html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.zip)
	}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	val(version)
	val(html)
	val(zip)
output:
	path("*.zip")
	path("*.html")
	path("*.xml")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
ln -s ${version} ./
ln -s ${html} ./
ln -s ${zip} ./
"""
}

runOnComplete(workflow);

