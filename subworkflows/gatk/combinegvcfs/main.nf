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
include { COMBINEGVCFS as HC_COMBINE1 } from '../../../modules/gatk/combinegvcfs'
include { COMBINEGVCFS as HC_COMBINE2 } from '../../../modules/gatk/combinegvcfs'
include { makeKey                     } from '../../../modules/utils/functions.nf'


List makeSQRT(def bed,def vcf_files) {
	def L = vcf_files.sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	if(n<25) n=25;
	def returnList = [];
	def currList = [];
	int i=0;
	for(;;) {
		if(i<L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) returnList.add([bed,currList]);
			if(i==L.size()) break;
			currList=[];
			}
		i++;
		}
	return returnList;
	}


workflow COMBINE_GVCFS {
take:
    meta
    fasta
    fai
    dict
    gvcfs_bed // [meta, gvcfgz, gvcfgz_tbi, bed]

main:
    versions= Channel.empty()
    ch1 = gvcfs_bed.map{meta,gvcf,tbi,bed ->
			[
			bed.toRealPath(),
			[gvcf,tbi]
			]}
		.groupTuple()
		.flatMap{bed,vcf_files->makeSQRT(bed,vcf_files)}
		.map{bed,vcf_files->[bed,vcf_files.flatten()]}
		.map{bed,vcf_files->[[id:makeKey([bed,vcf_files])],vcf_files,bed]}
	
	HC_COMBINE1(
		fasta,
		fai,
		dict,
		ch1
		)
	versions = versions.mix(HC_COMBINE1.out.versions)


	ch2 = HC_COMBINE1.out.gvcf
		.map{meta,gvcf,tbi,bed->[
			bed.toRealPath(),
			[gvcf,tbi]
		]}
		.groupTuple()
		.branch {bed,vcf_files->
			need_combine2:vcf_files.size()>1
			other:true /* only one gvcf for one bed, no need to combine a second level */
			}
	
	/** reformat channel for Combine2  */
	need_combine2 = ch2.need_combine2
		.map{bed,vcf_files->[bed,vcf_files.flatten().sort()]}
		.map{bed,vcf_files->[
			[id:makeKey([bed,vcf_files])],
			vcf_files,
			bed
			]}
	
	HC_COMBINE2(
		fasta,
		fai,
		dict,
		need_combine2
		)
	versions = versions.mix(HC_COMBINE2.out.versions)
	
	combined = ch2.other
		.map{bed,vcf_files->[bed,vcf_files.flatten()]}
		.map{bed,vcf_files->[
			vcf_files.find{it.name.endsWith(".vcf.gz")},
			vcf_files.find{it.name.endsWith(".vcf.gz.tbi")},
			bed
			]}
		.map{vcf,tbi,bed->[
			[id:makeKey(bed)],
			vcf,
			tbi,
			bed
			]}
		.mix(HC_COMBINE2.out.gvcf)

emit:
    versions
	gvcf = combined

}
