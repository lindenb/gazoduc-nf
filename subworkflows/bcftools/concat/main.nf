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
include {BCFTOOLS_CONCAT as CONCAT0 } from '../../../modules/bcftools/concat2'
include {BCFTOOLS_CONCAT as CONCAT1 } from '../../../modules/bcftools/concat2'
include {BCFTOOLS_CONCAT as CONCAT2 } from '../../../modules/bcftools/concat2'
include { makeKey                   } from '../../../modules/utils/functions.nf'
include { isBlank                   } from '../../../modules/utils/functions.nf'
include { verify                    } from '../../../modules/utils/functions.nf'



List makeSQRT(def meta,def vcf_files) {
	def L = vcf_files.sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	def returnList = [];
	def currList = [];
	int i=0;
	// save current group_id
	for(;;) {
		if(i < L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) {
				def h = [:];
				h =  h.plus(group_id: meta.id);
				h =  h.plus(id: meta.id+"_g"+i);
				returnList.add([h, currList]);
				}
			if(i==L.size()) break;
			currList=[];
			}
		i++;
		}
	return returnList;
	}

workflow BCFTOOLS_CONCAT {
take:
	metadata
	bed
	vcfs /* IMPORTANT : tuple [meta,vcf,idx] all VCF will be grouped by ID, so if you want to group everything, all the vcf will be grouped to one  */
main:
	if(metadata.min_group_size==null) {
		log.warn("undefined BCFTOOLS_CONCAT metadata.min_group_size");
		}
	def min_group_size = (metadata.min_group_size?:100)
	versions = Channel.empty()
	multiqc = Channel.empty()
	vcf_out = Channel.empty()


	/** group by id and count number of vcf per group */
	group1_ch = vcfs.map{meta,vcf,tbi->[meta.id,1]}
		.groupTuple()
		.map{group_id,array->[[id:group_id], (array.size() <= min_group_size) ]}

	ch1 = vcfs.combine(group1_ch)
		.filter{meta1,vcf,tbi,meta2,flag-> meta1.id == meta2.id}
		.branch {meta1,vcf,tbi,meta2,flag->
			level0  : flag == true //no need to split by sqrt
			level12 : true
			}

	CONCAT0(
		bed,
		ch1.level0.map{meta,vcf,tbi,_meta2,_flag->[meta.id,meta,[vcf,tbi]]}
			.groupTuple()
			.map{_group_id,metas,files->[ metas[0],files.flatten().sort()] }
		)
	versions = vcf_out.mix(CONCAT0.out.versions)
	vcf_out = vcf_out.mix(CONCAT0.out.vcf)
	

	ch2 = ch1.level12
		.map{meta,vcf,idx,_meta2,_flag->[[id:meta.id],[vcf,idx]]}
		.groupTuple()
		.flatMap{meta,vcf_files->makeSQRT(meta,vcf_files)}
		.map{meta,vcf_files->[meta, vcf_files.flatten().sort()]}
	

	CONCAT1(bed,ch2)
	versions = vcf_out.mix(CONCAT1.out.versions)
	vcf_out = vcf_out.mix(CONCAT1.out.vcf)

	
	ch3 = CONCAT1.out.vcf.view()
		.map{meta,vcf,tbi->[ meta.plus(id: meta.group_id)  , [ vcf,tbi]  ]}
		.groupTuple()
		.view()
		//.map{meta,vcf_files->[meta, vcf_files.flatten().sort()]}
		//.view()		
	/*
	CONCAT2(
		bed,
		ch2.map{meta,vcffiles->[
				meta,
				vcffiles.flatten().sort()
				]},
		)
	versions = versions.mix(CONCAT2.out.versions)
	*/
	
emit:	
	versions
	multiqc
	vcf = vcf_out //meta,vcf,tbi
}
