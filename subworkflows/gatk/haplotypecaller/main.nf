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

include { HAPLOTYPECALLER as HAPCALLER         }  from '../../../modules/gatk/hapcaller1'
include { COMBINEGVCFS as HC_COMBINE0          }  from '../../../modules/gatk/combinegvcfs'
include { COMBINEGVCFS as HC_COMBINE1          }  from '../../../modules/gatk/combinegvcfs'
include { GENOTYPEGVCFS                        }  from '../../../modules/gatk/genotypegvcfs'


List makeSQRT(def L1) {
	def key = L1.get(0);
	def L = L1.get(1).sort();
	int n = (int)Math.ceil(Math.sqrt(L.size()));
	if(n<25) n=25;
	def returnList = [];
	def currList = [];
	int i=0;
	for(;;) {
		if(i<L.size()) currList.add(L.get(i));
		if(i==L.size() || currList.size()==n) {
			if(!currList.isEmpty()) returnList.add([key,currList]);
			if(i==L.size()) break;
			currList=[];
			}
		i++;
		}
	return returnList;
	}


workflow HAPLOTYPECALLER {
take:
    meta
    fasta
    fai
    dict
    beds // meta,bed
    bams // meta,bam,bai
main:
    versions = Channel.empty()

    HAPCALLER(
        fasta,
        fai,
        dict,
        bams.combine(beds.map{it[1]})
        )
    versions  = versions.mix(HAPCALLER.out.versions)

    ch1 = HAPCALLER.out.gvcf
        .map{[[it[3].toRealPath()],[it[1],it[2]]]}
        .groupTuple()

    ch2 = ch1.branch{v->
        //nocombine: v[1].size()==1
        combine1 : v[1].size() < 100
        combine2 : true
        }

  to_genotype = Channel.empty()
  /***************************************************
   *
   * NO COMBINE NEEDED
   *
   */
//   to_genotype = to_genotype.mix(ch2.nocombine.map{[[id:it[0].toString().md5().substring(0,8)],it[1].flatten()]})
  
  /***************************************************
   *
   * LEVEL1
   *
   */    
  HC_COMBINE0(
	fasta,
	fai,
	dict,
	ch2.combine1.map{[[id:it[0].toString().md5().substring(0,7)],it[1].flatten(),it[0]]}
	)
  versions  = versions.mix(HC_COMBINE0.out.versions)
  to_genotype= to_genotype.mix(HC_COMBINE0.out.gvcf)
  /***************************************************
   *
   * LEVEL2
   *
   */    
   /*level2a_ch = ch2.combine2
            .flatMap{makeSQRT(it)}
			.map{[[id:it[0].name],it[1].flatten(),it[0]]}
    
    HC_COMBINE1(fasta,fai,dict, level2a_ch )
    versions  = versions.mix(HC_COMBINE1.out.versions)

*/
    to_genotype
    GENOTYPEGVCFS(
        fasta,
        fai,
        dict,
        to_genotype
        )
    versions  = versions.mix(GENOTYPEGVCFS.out.versions)
    
emit:
    versions
}
