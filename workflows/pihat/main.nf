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


include {dumpParams;runOnComplete            } from '../../modules/utils/functions.nf'
include {makeKey                             } from '../../modules/utils/functions.nf'
include {isBlank                             } from '../../modules/utils/functions.nf'
include {PIHAT as PIHAT01                    } from '../../subworkflows/pihat'
include {assertKeyExistsAndNotEmpty          } from '../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty            } from '../../modules/utils/functions.nf'
include {PREPARE_ONE_REFERENCE               } from '../../subworkflows/samtools/prepare.one.ref'
include {MULTIQC                             } from '../../subworkflows/multiqc'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()
	def workflow_meta = [
		id: "pihat"
		]

	PREPARE_ONE_REFERENCE(
		workflow_meta.plus([skip_scatter:true]),
		Channel.fromPath(params.fasta).map{[[id:it.baseName],it]}
		)
    versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/* load gnomad */
	if(params.gnomad==null) {
		gnomad_ch = Channel.empty()
		}
	else if(params.gnomad.endsWith(".vcf.gz")) {
		gnomad_ch = Channel.of(params.gnomad)
			.map{[[id:"gnomad"],file(it),file(it+".tbi")]}
		}
	else if(params.gnomad.endsWith(".list")) {
		gnomad_ch = Channel.fromPath(params.gnomad)
			.splitText()
			.map{it.trim()}
			.map{[[id:"gnomad"],file(it),file(it+".tbi")]}
		}
	else {
		gnomad_ch = Channel.fromPath(params.gnomad)
			.splitCsv(header:true,sep:',')
			.map{assertKeyExistsAndNotEmpty(it,"vcf")}
			.map{
				if(isBlank(it.tbi)) return it.plus("tbi":it.vcf+".tbi");
				return it;
				}
			.map{[[id:"gnomad"],file(it.vcf),file(it.tbi)]}
		}

	/* load 1000G */
	if(params.onekg==null) {
		onekg_ch = Channel.empty()
		}
	else if(params.onekg.endsWith(".vcf.gz")) {
		onekg_ch = Channel.of(params.onekg)
			.map{[[id:"1kg"],file(it),file(it+".tbi")]}
		}
	else if(params.onekg.endsWith(".list")) {
		onekg_ch = Channel.fromPath(params.onekg)
			.splitText()
			.map{it.trim()}
			.map{[[id:"1kg"],file(it),file(it+".tbi")]}
		}
	else {
		onekg_ch = Channel.fromPath(params.onekg)
			.splitCsv(header:true,sep:',')
			.map{assertKeyExistsAndNotEmpty(it,"vcf")}
			.map{
				if(isBlank(it.tbi)) return it.plus("tbi":it.vcf+".tbi");
				return it;
				}
			.map{[[id:"1kg"],file(it.vcf),file(it.tbi)]}
		}

	/* load user's VCF */
	if(params.vcf==null) {
		throw new IllegalArgumentException("undefined option --vcf.")
		}
	else if(params.vcf.endsWith(".vcf.gz")) {
		uservcf_ch = Channel.of(params.vcf)
			.map{[[id:makeKey(it)],file(it),file(it+".tbi")]}
		}
	else if(params.vcf.endsWith(".list")) {
		uservcf_ch = Channel.fromPath(params.vcf)
			.splitText()
			.map{it.trim()}
			.map{[[id:makeKey(it)],file(it),file(it+".tbi")]}
		}
	else {
		uservcf_ch = Channel.fromPath(params.vcf)
			.splitCsv(header:true,sep:',')
			.map{assertKeyExistsAndNotEmpty(it,"vcf")}
			.map{
				if(isBlank(it.tbi)) return it.plus("tbi":it.vcf+".tbi");
				return it;
				}
			.map{[[id:makeKey(it)],file(it.vcf),file(it.tbi)]}
		}

	
	if(params.sample2pop==null) {
		sample2pop_ch = Channel.of([[id:"nosample2pop"],[]])
		}
	else
		{
		sample2pop_ch = Channel.of([[id:"sample2pop"],file(params.sample2pop)])
		}

	if(params.remove_samples==null) {
		exclude_samples_ch = Channel.of([[id:"no_x_samples"],[]])
		}
	else
		{
		exclude_samples_ch = Channel.of([[id:"exclude_samples"],file(params.remove_samples)])
		}
	if(params.exclude_bed==null) {
		exclude_bed_ch = Channel.of([[id:"no_x_bed"],[]])
		}
	else
		{
		exclude_bed_ch = Channel.of([[id:"exclude_bed"],file(params.exclude_bed)])
		}

	pihat_ch = PIHAT01(
		workflow_meta.plus(level:1),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		sample2pop_ch,
		exclude_samples_ch,
		exclude_bed_ch,
		onekg_ch,
		gnomad_ch,
		uservcf_ch
		)
	multiqc = multiqc.mix(PIHAT01.out.multiqc)
	
	
	MULTIQC(
		workflow_meta.plus("id":"pihat"),
		sample2pop_ch,
		versions,
		[[id:"no_mqc_config"],[]],
		multiqc
		)
	}

runOnComplete(workflow);
