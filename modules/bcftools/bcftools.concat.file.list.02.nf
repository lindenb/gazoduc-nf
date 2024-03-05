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
include {moduleLoad} from '../../modules/utils/functions.nf'


process BCFTOOL_CONCAT_FILE_LIST {
tag "${vcfs.name}"
afterScript "rm -rf TMP"
cpus 1
input:
        path(vcfs)
	path(bed) // or NO_FILE
output:
        tuple path("*.{bcf,vcf.gz}"),path("*.{csi,tbi}"),emit:output
script:
	def use_vcf=false
	def args = "--no-version --allow-overlaps --remove-duplicates "
"""	
	hostname 1>&2
	${moduleLoad("bcftools")}
	mkdir -p TMP

	bcftools concat --threads ${task.cpus} ${bed.name.equals("NO_FILE")?"":"--regions-file \"${bed}\""} \\
		${args} \\
		-O ${use_vcf?"z":"b"} \\
		-o "TMP/${params.prefix?:""}concatenated.${use_vcf?"vcf.gz":"bcf"}" --file-list "${vcfs}"


	bcftools index --force ${use_vcf?"--tbi":""} --threads ${task.cpus}  "TMP/${params.prefix?:""}concatenated.${use_vcf?"vcf.gz":"bcf"}"

	mv TMP/* ./
"""
}
