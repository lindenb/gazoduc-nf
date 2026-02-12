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
include {BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_LEVEL1    } from '../../../modules/bcftools/concat3'
include {BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_LEVEL2    } from '../../../modules/bcftools/concat3'
include { makeKey                                     } from '../../../modules/utils/functions.nf'

List makeSQRT(def meta,def vcf_files) {
        int n = (int)Math.ceil(Math.sqrt(vcf_files.size()));
        if(n<25) n=25;
        def returnList = [];
        def currList = [];
        int i=0;
        for(;;) {
                if(i < vcf_files.size()) currList.add(vcf_files.get(i));
                if(i == vcf_files.size() || currList.size()==n) {
                        if(!currList.isEmpty()) returnList.add([meta,currList.flatten().sort()]);
                        if(i==vcf_files.size()) break;
                        currList=[];
                        }
                i++;
                }
        return returnList;
        }



/******
 * Usage: VCF file will be grouped by meta.id
 */
workflow BCFTOOLS_CONCAT_ALL {
take:
	metadata
	vcfs //meta,vcf,tbi
main:
	versions = Channel.empty()
	multiqc  = Channel.empty()
	
	leve1_ch = vcfs
		.map{meta,vcf,tbi->[meta,[vcf,tbi]]}
		.groupTuple()
		.map{meta,files->[meta,files.sort{f1,f2->f1[0].name<=>f2[0].name}]}
		.flatMap(row->makeSQRT(row))
		.map{meta,array->[meta.plus(save_id:meta.id,id:makeKey(array)),array]} //save id, create new id to avoid name collision in LEVEL2
        
	BCFTOOLS_CONCAT_LEVEL1(leve1_ch)
	versions = versions.mix(BCFTOOLS_CONCAT_LEVEL1.out.versions)
	
	ch2 = BCFTOOLS_CONCAT_LEVEL1.out.vcf
		.map{meta,vcf,tbi->[meta.plus(id:meta.save_id),[vcf,tbi]]}
		.groupTuple()
	
	out_vcf = Channel.empty()

	ch3 = ch2.branch{meta,array->
		it_is_done : array.size()==1
		need_level2 : true
		}
	// add to vcf, those who just need level1
	out_vcf = out_vcf.mix( ch3.it_is_done.map{meta,array->[meta.minus(save_id:meta.save_id),array[0][0],array[0][1]]})

	//apply level 2
	BCFTOOLS_CONCAT_LEVEL2(ch3.need_level2.map{meta,array->[meta.minus(save_id:meta.save_id), array.flatten().sort()]})
	versions = versions.mix(BCFTOOLS_CONCAT_LEVEL2.out.versions)
	out_vcf = out_vcf.mix( BCFTOOLS_CONCAT_LEVEL2.out.vcf)

	
emit:
	vcf = out_vcf
	versions
	multiqc
}

