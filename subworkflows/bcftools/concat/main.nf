/*

Copyright (c) 2025 Pierre Lindenbaum

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
include {BCFTOOLS_CONCAT as CONCAT1 } from '../../../modules/bcftools/concat'
include {BCFTOOLS_CONCAT as CONCAT2 } from '../../../modules/bcftools/concat'
include {makeKey                    } from '../../../modules/utils/functions.nf'


List makeSQRT(def group_id,def vcf_files) {
	def L = vcf_files.sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	if(n<25) n=25;
	def returnList = [];
	def currList = [];
	int i=0;
	for(;;) {
		if(i<L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) returnList.add([group_id,currList]);
			if(i==L.size()) break;
			currList=[];
			}
		i++;
		}
	return returnList;
	}

workflow BCFTOOLS_CONCAT {
take:
	meta
	bed //meta,bed
	vcfs // tuple [meta,vcf,idx]
main:

	versions = Channel.empty()
	
	ch1 = vcfs
		.map{meta,vcf,idx->[meta.id,[vcf,idx]]}
		.groupTuple()
		.flatMap{group_id,vcf_files->makeSQRT(group_id,vcf_files)}
		.map{[[id:it[0],old_id:it[0]], it[1].flatten().sort()]}

	CONCAT1(ch1,bed)
	versions = versions.mix(CONCAT1.out.versions)

	vcf_out = Channel.empty()

	ch2 = CONCAT1.out.vcf
		.map{meta,vcf,tbi->[meta,(vcf instanceof List?vcf:[vcf]),(tbi instanceof List?tbi:[tbi]) ]}
		.map{[it[0].id,it[1].plus(it[2])]}
		.groupTuple()
		.branch{
			done: it[1].size()==1 // ONE array containing [vcf,idx]
			need_concat2 : true
			}

	


	vcf_out= ch2.done.map{[
		[id:it[0]],
		it[1][0].find{f->f.name.endsWith(".bcf") || f.name.endsWith(".vcf.gz")},
		it[1][0].find{f->f.name.endsWith(".csi") || f.name.endsWith(".tbi")}
		]}


	CONCAT2(
		ch2.need_concat2
			.map{[
				[id:it[0]],
				it[1].flatten().sort()
				]},
		bed
		)
	versions = versions.mix(CONCAT2.out.versions)

	vcf_out = vcf_out.mix(
		CONCAT2.out.vcf
			.map{meta,vcf,tbi->[meta,(vcf instanceof List?vcf:[vcf]),(tbi instanceof List?tbi:[tbi]) ]}
			.flatMap{meta,vcfs,tbis->{
				def L=[];
				def X1 = vcfs.sort();
				def X2 = tbis.sort();
				for(int i=0;i< X1.size();i++) {
					L.add([meta,X1[i],X2[i]]);
					}
				return L;
				}}
		)
emit:	
	versions
	vcf = vcf_out
}