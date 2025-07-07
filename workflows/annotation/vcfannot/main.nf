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
nextflow.enable.dsl=2

include {runOnComplete} from '../../../modules/utils/functions.nf'
//include {ANNOTATE_VCF_01} from '../../../subworkflows/annotation/annotation.vcf.01.nf'
//include {MERGE_VERSION} from '../../../modules/version/version.merge.02.nf'
//include {JVARKIT_VCF_TO_INTERVALS_01} from '../../../subworkflows/jvarkit/vcf2intervals/main.nf'
//include {COLLECT_TO_FILE_01} from '../../../modules/utils/collect2file.01.nf'
//include {BCFTOOLS_CONCAT_01} from '../../../subworkflows/bcftools/bcftools.concat.01.nf'
//include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed/main.nf'
//include {VCF_TO_BED} from '../../../modules/bcftools/vcf2bed/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GTF   } from '../../../modules/gtf/download/main.nf'
include {DOWNLOAD_GTF_OR_GFF3 as DOWNLOAD_GFF3  } from '../../../modules/gtf/download/main.nf'
//include {ANNOTATE                               } from '../../../subworkflows/annotation/annotsnv/main.nf'
include {ANNOTATE                                 } from '../../../subworkflows/annotation/annotsnv/main.nf'
include {SCATTER_TO_BED                           } from '../../../subworkflows/gatk/scatterintervals2bed/main.nf'
include {SPLIT_N_VARIANTS                         } from '../../../modules/jvarkit/splitnvariants/main.nf'
include {BCFTOOL_CONCAT                           } from '../../../modules/bcftools/concat/main.nf'

workflow {
	versions = Channel.empty()
	def ref_hash = [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName
		]
	def fasta =    [ ref_hash, file(params.fasta) ]
	def fai =      [ ref_hash, file(params.fai) ]
	def dict =     [ ref_hash, file(params.dict) ]
	def pedigree = [ ref_hash, (params.pedigree ? file(params.pedigree) : []) ]
	def bed =      [ ref_hash, (params.bed ? file(params.bed) : []) ]
	
	


	vcfs = Channel.fromPath(params.samplesheet)
		.splitCsv(header:true,sep:",")
		.map{
			if(!it.vcf) throw new IllegalArgumentException("vcf missing in ${params.samplesheet}")
			if(!(it.vcf.endsWith(".vcf.gz") || it.vcf.endsWith(".bcf"))) throw new IllegalArgumentException("bad vcf suffix ${it.vcf} in ${params.samplesheet}")
			return it;
			}
		.map{
			if(it.containsKey("index")) return it;
			def suffix = it.vcf.endsWith(".vcf.gz")?".tbi":".csi"
			return it.plus(index:it.vcf+suffix);
			}
		.map{
			if(it.containsKey("id")) return it;
			return it.plus(id:file(it.vcf).baseName);
			}
		.map {[[id:it.id],file(it.vcf),file(it.index)]}
		.view()


	if(!bed[1]) {
		SCATTER_TO_BED([id:"annot"],fasta,fai,dict)
		bed = SCATTER_TO_BED.out.bed
		versions= versions.mix(SCATTER_TO_BED.out.versions)
	} 

	SPLIT_N_VARIANTS(
		vcfs.combine(bed).map{[it[0],it[1],it[2],it[4]]}
		)
	versions= versions.mix(SPLIT_N_VARIANTS.out.versions)

	ch1 = SPLIT_N_VARIANTS.out.vcf.
		flatMap{
			def L=[]
			for(X in it[1]) {
				L.add([it[0],X])
			}
		return L;
		}

	ch2= SPLIT_N_VARIANTS.out.tbi
		.flatMap{T->{
			def L=[]
			for(X in T[1]) {
				L.add([T[0],X])
				}
			return L;
			}}
	
	vcf = ch1.combine(ch2)
		.filter{it[0].equals(it[2])}
		.filter{(it[1].toRealPath()+".tbi").equals(it[3].toRealPath())}
		.map{[it[0],it[1],it[3]]}
		.view()


	DOWNLOAD_GTF(fasta,fai,dict)
	DOWNLOAD_GFF3(fasta,fai,dict)
	
	ANNOTATE(
		[id:"annot"],
		fasta,
		fai,
		dict,
		DOWNLOAD_GTF.out.output,
		DOWNLOAD_GFF3.out.output,		
		pedigree,
		vcf
		)

	BCFTOOL_CONCAT(
		ANNOTATE.out.vcf
			.map{[ it[0], [it[1],it[2]] ]}
			.groupTuple()
			.map{[it[0],it[1].flatten()]} , 
		[[id:"nobed"],[]]
		)

	versions= versions.mix(ANNOTATE.out.versions)
	}

workflow ANNOTATE_VCF {
	take:
		meta
		fasta
		fai
		dict
		vcf
		bed
		pedigree
	main:
		version_ch = Channel.empty()
		
		

		/*
		ann_ch = ANNOTATE(meta,
			fasta,
			fai,
			dict,
			bed,
			SPLIT_N_VARIANTS()
			)
		version_ch = version_ch.mix(ann_ch.versions)
		*/

	/*
		tofile_ch = COLLECT_TO_FILE_01([suffix:".list"], ann_ch.output.map{T->T.annot_vcf}.collect())
		version_ch = version_ch.mix(tofile_ch.versions)

		concat_ch = BCFTOOLS_CONCAT_01([:], tofile_ch.output, file("NO_FILE") )
		version_ch = version_ch.mix(concat_ch.version)
	*/
	emit:
		versions
		//vcf = concat_ch.vcf
	}

runOnComplete(workflow);
