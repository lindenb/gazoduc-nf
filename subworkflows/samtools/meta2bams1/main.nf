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
include {assertKeyExistsAndNotEmpty          } from '../../../modules/utils/functions.nf'
include {isBlank                             } from '../../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES                    } from '../../../modules/samtools/samples'

workflow META_TO_BAMS {
take:
	meta
	fasta
	fai
	metas
main:
	versions = Channel.empty()

	ch1 = metas
		.map{assertKeyExistsAndNotEmpty(it,"bam")}
		.branch{
			no_sample: it.sample==null || it.sample.isEmpty() || it.sample.equals(".")
			with_sample: true
			}

	ch2 = ch1.no_sample.map{[file(it.bam).toRealPath(),it]}
	
	SAMTOOLS_SAMPLES(
		fasta
			.mix(fai)
			.map{meta,f->f}
			.collect()
			.map{[[id:"ref"],it.sort()]},
		ch2
			.map{bam,meta->bam}
			.collect()
			.map{[[id:"bams"],it.sort()]}
		)
	versions = versions.mix(SAMTOOLS_SAMPLES.out.versions)

	
	ch3 = SAMTOOLS_SAMPLES.out.samplesheet
		.splitCsv(sep:'\t',header:false)
		.map{if(isBlank(it[0])) throw new IllegalArgumentException("No sample name for ${it}"); return it;}
		.map{if(it.size()< 2 || isBlank(it[2])) throw new IllegalArgumentException("Not same reference for ${it}"); return it;}
		.map{row->[file(row[1]).toRealPath().toString(),[sample:row[0]]]}
	
	
	fixed_sample_ch= ch2
		.map{bam,meta->[bam.toRealPath().toString(),meta]}
		.join(ch3, failOnDuplicate:true, failOnMismatch:true)
		.map{bam,meta1,meta2->meta2.plus(meta1)}
	
	bams_out = ch1.with_sample
		.mix(fixed_sample_ch)
        .map{if(!(it.bam.endsWith(".cram") || it.bam.endsWith(".bam"))) throw new IllegalArgumentException("${it}.bam should end with bam or cram"); return it;}
		.map{it.bai?it: (it.bam.endsWith(".bam") ? it.plus(["bai":it.bam+".bai"]):  it.plus(["bai":it.bam+".crai"]))}
		.map{it.id!=null?it:it.plus(id:it.sample)}
		.map{assertKeyExistsAndNotEmpty(it,"id")}
		.map{[
            it.findAll{k,v->!k.matches("(bam|bai|fasta|fai|dict)") && !isBlank(v)},
            file(it.bam),
            file(it.bai)
            ]};

	if(meta.with_test_unique_id==null || meta.with_test_unique_id==true) {
     	bams_out.map{it->it[0].id}.unique().count().
     		combine(bams_out.map{it->it[0].id}.unique().count())
     		.filter{c1,c2->c1!=c2}
     		.map{throw new IllegalArgumentException("input contains same duplicate samples");}
     	}

emit:
	versions
	bams = bams_out

}
