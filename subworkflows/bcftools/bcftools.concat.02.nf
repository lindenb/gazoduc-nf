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

include {BCFTOOL_CONCAT_FILE_LIST as CONCAT1} from '../../modules/bcftools/bcftools.concat.file.list.02.nf' addParams(prefix:"tmp.")
include {BCFTOOL_CONCAT_FILE_LIST as CONCAT2} from '../../modules/bcftools/bcftools.concat.file.list.02.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {SQRT_FILE} from '../../modules/utils/sqrt.nf'

workflow BCFTOOLS_CONCAT {
take:
	meta
	vcfs /* stream of hash map with key 'vcf' */
	bed /* or NO_FILE */
main:
	if(!meta.containsKey("method")) throw new IllegalArgumentException("concat.method is Missing. use method='all' for default.");

	if(meta.method.equals("all")) {
		vcflist1 = vcfs.map{it.vcf}.collectFile(name: 'to_concat.list', newLine: true,sort:'hash')
		d1_ch = SQRT_FILE(min_file_split:meta.min_file_split?:20, suffix:".list", vcflist1)
		d2_ch = d1_ch.output.splitText().map{file(it.trim())}

        	d3_ch = CONCAT1(d2_ch,bed)

		vcflist2 = d3_ch.output.map{T->T[0].toString()}.collectFile(name: 'to_concat.list', newLine: true,sort:'hash')

		d4_ch = CONCAT2(vcflist2 ,bed)

				
		out1 = d4_ch.output.map{T->[vcf:T[0][0],index:T[0][1]]}

		out2 =  d4_ch.output.map{T->T[0].toString()}.collectFile(name: 'vcfs.list', newLine: true,sort:'hash')
		}
	else
		{
		throw new IllegalArgumentException("concat.method unknown ${meta.method}");
		}

emit:
	vcfs = out1 /* stream of hash [vcf:,index:]
	vcf_list = out2 /* vcfs as list */
}



