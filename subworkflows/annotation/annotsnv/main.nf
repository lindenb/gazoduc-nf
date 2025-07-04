/*

Copyright (c) 2024 Pierre Lindenbaum

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

include {ALPHAMISSENSE} from '../alphamissense/main.nf'
include {BCFTOOLS_BCSQ} from  '../../../modules/bcftools/bcsq/main.nf'

/*
include {slurpJsonFile;moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {BCFTOOLS_CONCAT} from '../bcftools/bcftools.concat.02.nf'
include {STEP_FIRST} from './step.first.nf'
include {ANNOTATE_RMSK} from './step.rmsk.nf'
include {ANNOTATE_BCSQ} from './step.bcsq.nf'
include {ANNOTATE_CLINVAR} from './step.clinvar.nf'
include {ANNOTATE_VISTA} from './step.vista.nf'
include {ANNOTATE_SNPEFF} from './step.snpeff.nf'
include {ANNOTATE_SIMPLE_REPEATS} from './step.simple_repeats.nf'
include {ANNOTATE_REMAP} from './step.remap.nf'
include {ANNOTATE_GNOMADSV} from './step.gnomadsv.nf'
include {ANNOTATE_REGFEATURES} from './step.regfeatures.nf'
include {ANNOTATE_AVADA} from './step.avada.nf'
include {ANNOTATE_VEP} from './step.vep.nf'
include {ANNOTATE_GNOMAD_GENOME} from './step.gnomad.genome.nf'
include {ANNOTATE_CONTRAST} from './step.contrast.nf'
include {ANNOTATE_POLYX} from './step.polyx.nf'
include {ANNOTATE_ELSEWHERE} from './step.elsewhere.nf'
include {ANNOTATE_FILTERSO} from './step.vcffilterso.nf'
include {ANNOTATE_GENCC} from './step.gencc.nf'
include {ANNOTATE_GREENDB} from './step.greendb.nf'
include {ANNOTATE_BHFUCL} from './step.bhfucl.nf'
include {ANNOTATE_NORM} from './step.norm.nf'
include {ANNOTATE_SPLICEAI} from './step.spliceai.nf'
include {ANNOTATE_STRINGDB} from './step.stringdb.nf'
include {ANNOTATE_MONDO} from './step.mondo.nf'
include {ANNOTATE_CCRE} from './step.ccre.nf'
include {ANNOTATE_CADD} from './step.cadd.nf'
include {ANNOTATE_REGULOME} from './step.regulome.nf'
include {ANNOTATE_UORFDB} from './step.uorfdb.nf'
include {ANNOTATE_REVEL} from './step.revel.nf'
include {ANNOTATE_MNV} from './step.mnv.nf'
*/
/*

	step.gnomad.constraint.nf
	step.gtex.tmp.nf
	step.hmc.nf
	step.hpo.nf
	step.mcap.nf
	step.phastcons.nf
	step.trio-dnm2.nf
	step.uniprot.nf
*/
boolean hasFeature(hash,key) {
	def k2 = "with_"+key;
	if(!hash.containsKey(k2)) return false;
	return hash[k2].equals(true);
	}

workflow ANNOTATE {
	take:
		meta
		fasta
		fai
		dict
		gtf
		gff3
		pedigree
		vcf
	main:
		versions = Channel.empty()
		def NO_BED = [meta,[]];

		ALPHAMISSENSE(meta, fasta, fai, dict, NO_BED,vcf)
		versions = versions.mix(ALPHAMISSENSE.out.versions)
		vcf = ALPHAMISSENSE.out.vcf

		BCFTOOLS_BCSQ(fasta, fai, dict, gff3, vcf)
		versions = versions.mix(BCFTOOLS_BCSQ.out.versions)
		vcf = BCFTOOLS_BCSQ.out.vcf
		/*
		step_ch = STEP_FIRST(rows)
		step_ch = ANNOTATE_RMSK(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_VISTA(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_BCSQ(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_CLINVAR(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_SIMPLE_REPEATS(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_SNPEFF(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_FILTERSO(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_REMAP(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_GNOMAD_GENOME(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_GNOMADSV(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_CONTRAST(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
	
		step_ch = ANNOTATE_POLYX(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_AVADA(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_REGFEATURES(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)			
		
		step_ch = ANNOTATE_ELSEWHERE(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		

		step_ch = ANNOTATE_GENCC(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_GREENDB(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_BHFUCL(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
				
		step_ch = ANNOTATE_STRINGDB(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		

		step_ch = ANNOTATE_SPLICEAI(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_MONDO(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_CCRE(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_CADD(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_REGULOME(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_UORFDB(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_REVEL(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)


		step_ch = ANNOTATE_MNV(genomeId, bed, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		

		vcfs2 = step_ch.output.map{T->slurpJsonFile(T)}.map{T->[vcf:T.vcf,index:T.index]}
                

                concat_ch = BCFTOOLS_CONCAT([method:"all"], vcfs2 , file("NO_FILE"))
		

		version_ch = MERGE_VERSION("VCF annotation", version_ch.collect())
		*/
	emit:
		versions = versions_ch
	}


