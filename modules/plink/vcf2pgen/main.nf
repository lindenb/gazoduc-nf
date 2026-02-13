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
process PLINK2_VCF2PGEN {
label "process_single_high"
afterScript "rm -rf TMP"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta) /* or any ... */) //for PAR region
	tuple val(meta2),path(plink_ped)
	tuple val(meta ),path(vcf)
output:
	tuple val(meta),path("*"),emit:output
    path("versions.yml")
script:
	def k1 = k1_signature()
	def args = contig.matches("(chr)?[XY]")?" --update-sex ${plink_ped} --split-par ${meta1.ucsc_name?:"undefined"}":""
    def prefix = task.ext.prefix?:"${meta.id}"
"""
mkdir -p TMP

plink2 ${vcf.name.endsWith(".bcf")?"--bcf":"--vcf"} "${vcf}"  \\
	${args} \\
	--make-pgen erase-phase \\
	--threads ${task.cpus} \\
	--out "${prefix}"

touch versions.yml
"""
stub:
"""
touch versions.yml
"""
}
