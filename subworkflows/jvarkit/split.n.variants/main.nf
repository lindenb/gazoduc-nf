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

include {SPLIT_N_VARIANTS  as SPLIT_N_VAR    } from '../../../modules/jvarkit/splitnvariants/main.nf'
include { VCF_TO_BED                         } from '../../../modules/bcftools/vcf2bed'
include { BCFTOOLS_INDEX                     } from '../../../modules/bcftools/index'
include { LINUX_SPLIT                        } from '../../../modules/utils/linuxsplit'
include { SPLIT_BED_PER_CHROMOSOME           } from '../../../modules/bed/splitchr'
include { flatMapByIndex                     } from '../../../modules/utils/functions.nf'


workflow SPLIT_N_VARIANTS {
take:
    metadata
    input_bed //one or more bed files, e.g. output of bedcluster
    vcfs
main:
    versions = Channel.empty()
    multiqc = Channel.empty()
	

	vcfs.map{_meta,_vcf,_tbi->[_meta.id,_meta]}
		.groupTuple()
		.filter{_id,metas->metas.size()!=1}
		.map{_id,metas->
			throw new IllegalArgumentException("Multiple VCF with the same id e.g: ${metas}");
			return -1;
			}


	
	/** split input BED per chromosome */
	SPLIT_BED_PER_CHROMOSOME(input_bed)
	versions = versions.mix(SPLIT_BED_PER_CHROMOSOME.out.versions)

	contig_beds = SPLIT_BED_PER_CHROMOSOME.out.beds
		.map{meta,bed->[meta,(bed instanceof List?bed:[bed])]}
		.flatMap{_meta,beds->beds}
		.map{bed->[bed,bed]}//duplicate bed data
		.splitCsv(header:false,sep:'\t',limit:1) /* limit to one row per bed */
		.map{bed_row,bed->[bed_row[0]/* contig */,bed]}
		.unique()
		.map{contig,bed->[[id:bed.baseName,chrom:contig],bed]}
		

	VCF_TO_BED(vcfs)
	versions = versions.mix(VCF_TO_BED.out.versions)

	//vcfbed_ch = VCF_TO_BED.out.output.map{meta,bed,vcf,tbi->[meta.plus(id:makeKey(bed),src_id:bed),bed,vcf,tbi]}


	LINUX_SPLIT(VCF_TO_BED.out.output.map{meta,bed,_vcf,_tbi->[meta,bed]})
	versions = versions.mix(LINUX_SPLIT.out.versions)

	vcf_contig_ch = LINUX_SPLIT.out.output
		.map{meta,f->[meta,(f instanceof List?f:[f])]}
		.flatMap{flatMapByIndex(it,1)}
		.splitCsv(header:false,sep:'\t',limit:1)
		.map{meta,bed_row->[meta.plus(chrom:bed_row[0])]}
		.combine(VCF_TO_BED.out.output.map{meta,_bed,vcf,tbi->[meta,vcf,tbi]})
		.filter{meta1,meta2,_vcf,_tbi->meta1.id==meta2.id}
		.map{meta1,_meta2,vcf,tbi->[meta1,vcf,tbi]}

	

	ch1 = vcf_contig_ch
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

	BCFTOOLS_INDEX(
		SPLIT_N_VAR.out.vcf.flatMap{flatMapByIndex(it,1)}
		)
	versions= versions.mix(BCFTOOLS_INDEX.out.versions)
	
	BCFTOOLS_INDEX.out.vcf

	/* remove "meta.chrom' to keep original meta */ 
	vcf_out =BCFTOOLS_INDEX.out.vcf
		.map{meta,vcf,idx->[meta.findAll{k,v->!k.matches("(chrom)")},vcf,idx]}
	
emit:
    versions
	multiqc
    vcf = vcf_out
}
