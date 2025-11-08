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
include {BCFTOOLS_CONCAT as CONCAT1 } from '../../../modules/bcftools/concat2'
include {BCFTOOLS_CONCAT as CONCAT2 } from '../../../modules/bcftools/concat2'
include {makeKey                    } from '../../../modules/utils/functions.nf'
include {JOIN_VCF_TBI               } from '../../../subworkflows/join.vcf.tbi'


List makeSQRT(def min_group,def meta,def vcf_files) {
	def L = vcf_files.sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	if(n<min_group) n=min_group;
	def returnList = [];
	def currList = [];
	int i=0;
	for(;;) {
		if(i<L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) returnList.add([meta,currList]);
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
	def min_group = (meta.min_group_size?:25)
	versions = Channel.empty()


	ch1 = vcfs
		.map{meta,vcf,idx->[meta,[vcf,idx]]}
		.groupTuple()
		.flatMap{meta,vcf_files->makeSQRT(min_group,meta,vcf_files)}
		.map{meta,vcf_files->[meta, vcf_files.flatten().sort()]}
		

	CONCAT1(bed,ch1)
	versions = versions.mix(CONCAT1.out.versions)

	

	ch2 = CONCAT1.out.vcf
		.map{meta,vcf,tbi->[meta,(vcf instanceof List?vcf:[vcf]),(tbi instanceof List?tbi:[tbi]) ]}
		.map{meta,vcfs,tbis->[meta.id,vcfs.plus(tbis)]}
		.groupTuple()
		

	CONCAT2(
		bed,
		ch2.map{meta,vcffiles->[
				meta,
				vcffiles.flatten().sort()
				]},
		)
	versions = versions.mix(CONCAT2.out.versions)

	vcf_out = CONCAT2.out.vcf
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
		
	
emit:	
	versions
	vcf = vcf_out //meta,vcf,tbi
}
