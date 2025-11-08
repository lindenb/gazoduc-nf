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


include {WGSIM  as APPLY_WGSIM    } from '../../../modules/samtools/wgsim'



/**
 * Generate FASQ files from FASTA
 *
 */
workflow WGSIM {
	take:
		metadata
		fasta // [meta,fasta]
	main:
		versions = Channel.empty()
		
		if(metadata.random_seed==null) {
			log.warn("WGSIM: metadata.random_seed was not set.");
			}
		if(metadata.n_samples==null) {
			log.warn("WGSIM: metadata.n_samples was not set. Using default");
			}

		def rnd = new java.util.Random(((metadata.random_seed?:0L) as long));

	
	def n_samples= ((metadata.n_samples?:10) as int)
	ch1 = Channel.of(0..<n_samples)
			.map{n->[id:"S"+(n+1)]}
	
	ch1 = ch1
		.collate(3)
		.map{
			if(it.size()!=3) return it;
			def c = it[0].plus(["sex":(rnd.nextBoolean()?"male":"female")  ,"father":it[1].id, "mother":it[2].id])
			def f = it[1].plus(["sex":"male"  ,"father":"0",      "mother":"0"])
			def m = it[2].plus(["sex":"female","father":"0",      "mother":"0"])
			return [c,f,m];
			}
		.flatMap()
		.map{h->h.plus("status":(rnd.nextBoolean()?"case":"control"))}
		.map{h->h.plus("collection":"collection"+rnd.nextInt(3))}
		
		
	
	ch2 = ch1.combine(fasta)
		.map{meta1,meta2,fasta->[meta1.plus(reference_id:meta2.id),meta2,fasta]}
		.multiMap{meta1,meta2,fasta->
			fasta: [meta2,fasta]
			individual: meta1
			}
	
	APPLY_WGSIM(ch2.fasta,ch2.individual)
	versions = versions.mix(APPLY_WGSIM.out.versions)

	emit:
		versions
		fastq = APPLY_WGSIM.out.fastq
	}



