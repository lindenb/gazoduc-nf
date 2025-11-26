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

include {SPLIT_N_VARIANTS  as SPLIT_N_VAR    } from '../../../modules/jvarkit/splitnvariants/main.nf'
include { VCF_TO_BED                         } from '../../../modules/bcftools/vcf2bed'
include { SPLIT_BED_PER_CHROMOSOME           } from '../../../modules/bed/splitchr'

List toFlatMap(v) {
	def L1=[];
	def L2 = v[1];
	if(!(L2 instanceof List)) L2=[L2]
	for(f in L2) {
		L1.add([v[0],f])
		}
	return L1;
	}

workflow SPLIT_N_VARIANTS {
take:
    meta
    bed //one or more bed files, e.g. output of bedcluster
    vcfs
main:
    versions = Channel.empty()
    multiqc = Channel.empty()


	vcfs.map{_meta,_vcf,_tbi->[_meta.id,meta]}
		.groupTuple()
		.filter{id,metas->metas.size()!=1}
		.map{_id,metas->
			throw new IllegalArgumentException("Multiple VCF with the same id e.g: ${metas}");
			}

	SPLIT_BED_PER_CHROMOSOME(bed)
	versions = versions.mix(SPLIT_BED_PER_CHROMOSOME.out.versions)
	contig_beds = SPLIT_BED_PER_CHROMOSOME.out.beds
		.flatMap{meta,beds->beds}
		.map{bed->[bed,bed]}//duplicate bed data
		.splitCsv(header:false,sep:'\t',limit:1) /* limit to one row per bed */
		.map{bed_row,bed->[bed_row[0]/* contig */,bed]}
		.unique()
		.map{contig,bed->[[chrom:contig],bed]}

	VCF_TO_BED(vcfs)
	versions = versions.mix(VCF_TO_BED.out.versions)
	contig_vcfs = VCF_TO_BED.out.output
		.splitCsv(header:false,sep:'\t',limit:1) /* limit to one row per bed */)
		.map{meta,bed_row,vcf,tbi->[meta.plus(chrom:bed_row[1]),vcf,tbi]}
		.unique()

	ch1 = contig_vcfs
		.combine(contig_beds)
		.filter{meta1,_vcf,_tbi,meta2,_bed->meta1.chrom==meta2.chrom}
		.multiMap{meta1,vcf,tbi,meta2,bed->
			bed: [meta2,bed]
			vcf: [meta1,vcf,tbi]
			}

	SPLIT_N_VAR(
		ch1.bed,
		ch1.vcf
		)
	versions= versions.mix(SPLIT_N_VAR.out.versions)

	ch1 = SPLIT_N_VAR.out.vcf.
		flatMap{toFlatMap(it)}

	ch2= SPLIT_N_VAR.out.tbi
		.flatMap{toFlatMap(it)}
	
	vcf = ch1.combine(ch2)
		.filter{it[0].equals(it[2])}
		.filter{(it[1].toRealPath()+".tbi").equals(it[3].toRealPath())}
		.map{[it[0],it[1],it[3]]}
emit:
    versions
	multiqc
    vcf
}
