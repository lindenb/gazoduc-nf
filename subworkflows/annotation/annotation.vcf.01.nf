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
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
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
workflow ANNOTATE_VCF_01 {
	take:
		genomeId
		rows /* vcf, idx, bed, interval */
	main:
		version_ch = Channel.empty()
		step_ch = STEP_FIRST(rows)
		step_ch = ANNOTATE_RMSK(genomeId, step_ch.output)
		step_ch = ANNOTATE_VISTA(genomeId, step_ch.output)
		step_ch = ANNOTATE_BCSQ(genomeId, step_ch.output)
		step_ch = ANNOTATE_CLINVAR(genomeId, step_ch.output)
		step_ch = ANNOTATE_SIMPLE_REPEATS(genomeId, step_ch.output)
		step_ch = ANNOTATE_SNPEFF(genomeId, step_ch.output)
		step_ch = ANNOTATE_REMAP(genomeId, step_ch.output)
		step_ch = ANNOTATE_GNOMAD_GENOME(genomeId, step_ch.output)
		step_ch = ANNOTATE_GNOMADSV(genomeId, step_ch.output)

/*
		step_ch = ANNOTATE_REGFEATURES(genomeId, step_ch.output)
		step_ch = ANNOTATE_AVADA(genomeId, step_ch.output)
		step_ch = ANNOTATE_VEP(genomeId, step_ch.output)
*/
		
		rows=Channel.empty()

		version_ch = MERGE_VERSION("VCF annotation", version_ch.collect())
	emit:
		version= version_ch
		output = rows
	}


