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

include {moduleLoad;getKeyValue;hasFeature} from '../../../modules/utils/functions.nf'
include {WGSELECT_02} from '../../../subworkflows/wgselect/wgselect.02.nf'
include {VERSION_TO_HTML} from '../../../modules/version/version2html.nf'


params.reference=""
params.pedigree="NO_FILE"
params.vcf=""
params.disableFeatures="";
params.help=false
params.bed= "NO_FILE"

if(params.help) {
  log.info"""
## About

Cleanup variants from a VCF

## Author

Pierre Lindenbaum

## Options

  * --reference (fasta) The full path to the indexed fasta reference genome. It must be indexed with samtools faidx and with picard CreateSequenceDictionary or samtools dict. [REQUIRED]
  * --vcf <file> path to a indexed VCF or BCF file. If file ends with '.list' is a list of path to one VCF per contig [REQUIRED]
  * --pedigree <file> jvarkit formatted pedigree. phenotype MUST be case|control. Sex MUST be male|female|unknown
  * --bed <file> optional bed file to limit the analysis to the genes overlapping a  bed file.
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume wgselect.nf \\
        --publishDir output \\
        --prefix "analysis." \\
        --reference /path/to/reference.fasta \\
        --vcf /path/to/my.vcf.gz \\
        --pedigree /path/to/input.ped \
```

## Workflow

![workflow](./workflow.svg)

"""
exit 0
}

workflow {
		ch1 = WGSELECT_02(params, params.reference, file(params.vcf), file(params.pedigree), file(params.bed))
		html = VERSION_TO_HTML(params,ch1.version)
		PUBLISH(params, ch1.contig_vcfs, ch1.variants_list , ch1.version, html.html)
		}

process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
executor "local"
input:
	val(meta)
	path(vcfs)
	path(variants)
	path(xml)
	path(html)
output:

when:
        !params.getOrDefault("publishDir","").trim().isEmpty()
script:
	prefix = meta.prefix?:""
"""
mkdir "${prefix}archive"

for F in "${vcfs}" "${variants}" "${xml}" "${html}"
do
	ln -s "\${F}" "./${prefix}archive/" 
done

zip -r "${prefix}archive.zip" "${prefix}archive"
"""
}
