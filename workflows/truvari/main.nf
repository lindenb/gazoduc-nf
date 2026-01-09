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
nextflow.enable.dsl=2

include { PREPARE_ONE_REFERENCE       } from '../../subworkflows/samtools/prepare.one.ref'
include { READ_SAMPLESHEET            } from '../../subworkflows/nf/read_samplesheet'
include { assertKeyExistsAndNotEmpty  } from '../../modules/utils/functions.nf'
include { runOnComplete               } from '../../modules/utils/functions.nf'
include { removeCommonSuffixes        } from '../../modules/utils/functions.nf'
include { TRUVARI                     } from '../../subworkflows/truvari'
include { DOWNLOAD_GNOMAD_SV          } from '../../modules/gnomad_sv/download.vcf'
include { JVARKIT_VCFGNOMADSV         } from '../../modules/jvarkit/vcfgnomadsv'
include { BCFTOOLS_INDEX              } from '../../modules/bcftools/index'
include { BCFTOOLS_VIEW               } from '../../modules/bcftools/view'


workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()



	if(params.fasta==null) {
		throw new IllegalArgumentException("undefined --fasta");
		}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}

	def metadata = [id:"truvari"]
		

	PREPARE_ONE_REFERENCE(
		metadata,
		Channel.fromPath(params.fasta).map{f->[[id:f.baseName],f]}
		)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)
	
	READ_SAMPLESHEET(
		[arg_name:"samplesheet"],
		params.samplesheet
		)
	versions = versions.mix(READ_SAMPLESHEET.out.versions)
	vcf_ch = READ_SAMPLESHEET.out.samplesheet
		.map{row->assertKeyExistsAndNotEmpty(row,"vcf")}
		.map{row->file(row.vcf)}
		.map{f->[[id:removeCommonSuffixes(f.name)], f]}
	
	if(params.bed!=null) {
		BCFTOOLS_VIEW(
			[[id:"bed"],file(params.bed)],
			[[id:"nosample"],[]],
			vcf_ch.map{m,f->[m,f,[]]}
			)
		versions = versions.mix(BCFTOOLS_VIEW.out.versions)
		vcf_ch = BCFTOOLS_VIEW.out.vcf
		}


	DOWNLOAD_GNOMAD_SV( PREPARE_ONE_REFERENCE.out.dict )
    versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)
	
	JVARKIT_VCFGNOMADSV(
		DOWNLOAD_GNOMAD_SV.out.vcf,
		vcf_ch
		)
	versions = versions.mix(JVARKIT_VCFGNOMADSV.out.versions)
	vcf_ch = JVARKIT_VCFGNOMADSV.out.vcf


	BCFTOOLS_INDEX(JVARKIT_VCFGNOMADSV.out.vcf)
	versions = versions.mix(BCFTOOLS_INDEX.out.versions)
	vcf_ch = BCFTOOLS_INDEX.out.vcf

	TRUVARI(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		vcf_ch
		)
	versions = versions.mix(TRUVARI.out.versions)
	}

runOnComplete(workflow)
