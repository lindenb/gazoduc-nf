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
include { SEQKIT_SPLIT2   } from '../../../modules/seqkit/split2'
def restorePairs(def meta, def fastqs) {

	def L=[];
	int i=0;
	//System.err.println("restorePairs "+fastqs);

	def R1 = fastqs.findAll{F->F.name.contains("_R1_part_")}.sort();
	def R2 = fastqs.findAll{F->F.name.contains("_R2_part_")}.sort();

	if(R1.size()!=R2.size()) throw new IllegalArgumentException("restorePairs : R1.size != R2.size");

	if(R1.size()+R2.size()!=fastqs.size()) throw new IllegalArgumentException("restorePairs : R1.size+R2.size!=fastqs.size()");

	if(R1.isEmpty() && !fastqs.isEmpty()) throw new IllegalArgumentException("restore R1.size==0");

	for(i=0;i< R1.size();i++) {
		def hash = meta.plus([splid_id:"split_"+(i+1)])
		def f1 = R1[i];
		def f2 = R2[i];
		if(!f2.name.replaceAll("_R2_part_","_R1_part_").equals(f1.name)) throw new IllegalArgumentException("restorePairs faild: "+f1.name+" vs "+f2.name);
		L.add([hash, f1, f2]);
		}
	
	return L;
	}

def addSplitIndexForSingleEnd(def meta,def fastqs) {
	def L=[];
	int i=0;
	def R0 = fastqs.findAll{F->F.name.contains("_R0_part_")}.sort();
	if(R0.size() != fastqs.size()) throw new IllegalArgumentException("restore R0.size!=fastqs.size()");
	for(i=0;i< R0.size();i++) {
		def hash = meta.plus([splid_id:"split_"+(i+1)])
		L.add([hash, R0[i]]);
		}
	return L;
	}

/** split FASTQ using seqtk split2 */
workflow SEQKIT_SPLIT {
	take:
		meta
        fastqs // [meta, [R1,R2]] or [meta,R0] or [meta,R1,R2] or [meta,[R0]]
	main:		
		versions = Channel.empty()

		
  		fastqs = fastqs
                .map{
                    if(it.size()==2 && (it[1] instanceof Path)) {
                        return [it[0],it[1],[]];
                        }
					 if(it.size()==2 && (it[1] instanceof List) && it[1].size()==1) {
                        return [it[0],it[1][0],[]];
                        }
					if(it.size()==2 && (it[1] instanceof List) && it[1].size()==2) {
                        return [it[0],it[1][0], it[1][1]];
                        }
                    return it;
                    }
                .map {
                    if(it.size()!=3 || !(it[1] instanceof Path)) {
                        throw new IllegalArgumentException("workflow:SEQKIT_SPLIT: Bad input ${it}.");
                        }
                    return it;
                    }

        fastqs_out = fastqs

		if(meta.disable_split==null || meta.disable_split==false) {
			
			 

			SEQKIT_SPLIT2(fastqs)
			versions = versions.mix(SEQKIT_SPLIT2.out.versions)
			
			ch1 = SEQKIT_SPLIT2.out.fastqs
				.branch{_meta,fastqs->
					paired: fastqs.size()>1 &&  
							fastqs.findAll{F->F.name.contains("_R0_part_")}.size() ==0 && 
							fastqs.findAll{F->F.name.contains("_R1_part_")}.size() ==  fastqs.findAll{F->F.name.contains("_R2_part_")}.size()
					single: fastqs.findAll{F->F.name.contains("_R1_part_")}.size() == 0 && 
							fastqs.findAll{F->F.name.contains("_R2_part_")}.size() == 0 && 
							fastqs.findAll{F->F.name.contains("_R0_part_")}.size() ==  fastqs.size()
					others : true
				}

			ch1.others.map{throw new IllegalArgumentException("SEQTK SPLIT: boum ${it}.")};

			fastqs_out = Channel.empty()
			fastqs_out = fastqs_out.mix(
					ch1.single.flatMap{addSplitIndexForSingleEnd(it)}
					)

			fastqs_out = fastqs_out.mix(
					ch1.paired.flatMap{meta,fastqs->restorePairs(meta,fastqs)}
				)
		}

	emit:
		versions
		fastqs = fastqs_out
	}
