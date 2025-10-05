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
include {SEQKIT_SPLIT2 as SPLIT2            } from '../../../modules/seqkit/split2'

def restorePairs(def meta, def fastqs) {
	def L=[];
	def R1 = fastqs.findAll{F->F.name.endsWith("R1.fq.gz")}.sort()
	def R2 = fastqs.findAll{F->F.name.endsWith("R2.fq.gz")}.sort()
	if(R1.size()!=R2.size()) throw new IllegalArgumentException("restorePairs : R1.size != R2.size");
	if(R1.size()+R2.size()!=fastqs.size()) throw new IllegalArgumentException("restorePairs : R1.size+R2.size!=fastqs.size()");
	if(R1.isEmpty() && !fastqs.isEmpty()) throw new IllegalArgumentException("restore R1.size==0");
	for(int i=0;i< R1.size();i++) {
		def hash =  meta.plus([split_index:i+1]);
		def f1 = R1[i];
		def f2 = R2[i];
		if(!f2.name.replaceAll("R2.fq.gz").equals(f1.name)) throw new IllegalArgumentException("restorePairs : "+f1.name+" vs "+f2.name);
		L.add([hash, f1, f2]);
		}
	return L;
	}

/** split FASTQ using seqtk split2 */
workflow SEQTK_SPLIT {
	take:
		meta
		fastqs // meta, R1,optional_R2
	main:		
		versions = Channel.empty()
		
		ch1 = fastqs.branch{
			all: meta.sektq_split_args!=null && !meta.sektq_split_args.trim().isEmpty()
			seqtk: it[0].sektq_split_args!=null && !it[0].sektq_split_args.trim().isEmpty()
			other: true
			}


		out_fq  = ch1.other
		

		SPLIT2(ch1.seqtk.mix(ch1.all.map{meta,R1,R2->[meta.plus([sektq_split_args:meta.sektq_split_args]),R1,R2]}))
		versions = versions.mix(SPLIT2.out.versions)
		ch2 = SPLIT2.out.fastqs
			.view()
			.flatMap{meta,fastqs->restorePairs(meta,fastqs)}

		out_fq = out_fq.mix(ch2)
			.map{meta,R1,R2->[meta.findAll{k,v->!k.equals("sektq_split_args")},R1,R2]}
	emit:
		versions
		fastqs = out_fq
	}
