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
include {ANNOTATE_ALPHAMISSENSE} from './step.alphamissense.nf'
include {ANNOTATE_MONDO} from './step.mondo.nf'
include {ANNOTATE_CCRE} from './step.ccre.nf'

workflow ANNOTATE_VCF_01 {
	take:
		genomeId
		elsewhere_vfs
		rows /* vcf, idx, bed, interval */
	main:
		version_ch = Channel.empty()
		count_ch = Channel.empty()
		doc_ch = Channel.empty()
		
		step_ch = STEP_FIRST(rows)
		step_ch = ANNOTATE_RMSK(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_VISTA(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_BCSQ(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_CLINVAR(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_SIMPLE_REPEATS(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_SNPEFF(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		


		step_ch = ANNOTATE_FILTERSO(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_REMAP(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_GNOMAD_GENOME(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_GNOMADSV(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_CONTRAST(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
	
		step_ch = ANNOTATE_POLYX(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_AVADA(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_REGFEATURES(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)			
		
		step_ch = ANNOTATE_ELSEWHERE(genomeId, step_ch.output, elsewhere_vfs)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		

		step_ch = ANNOTATE_GENCC(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_GREENDB(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)		
		
		step_ch = ANNOTATE_BHFUCL(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		
		step_ch = ANNOTATE_STRINGDB(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		
		step_ch = ANNOTATE_ALPHAMISSENSE(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_SPLICEAI(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_MONDO(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)

		step_ch = ANNOTATE_CCRE(genomeId, step_ch.output)
		count_ch = count_ch.mix(step_ch.count)
		doc_ch = count_ch.mix(step_ch.doc)
		

		vcfs2 = step_ch.output.map{T->slurpJsonFile(T)}.map{T->[vcf:T.vcf,index:T.index]}
                

                concat_ch = BCFTOOLS_CONCAT([method:"all"], vcfs2 , file("NO_FILE"))
		

		version_ch = MERGE_VERSION("VCF annotation", version_ch.collect())
	emit:
		version= version_ch
		output = rows
	}


