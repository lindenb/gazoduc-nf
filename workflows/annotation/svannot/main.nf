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


include { validateParameters                          } from 'plugin/nf-schema'
include { paramsHelp                                  } from 'plugin/nf-schema'
include { paramsSummaryLog                            } from 'plugin/nf-schema'
include { samplesheetToList                           } from 'plugin/nf-schema'
include { runOnComplete;dumpParams                    } from '../../../modules/utils/functions.nf'
include { PREPARE_ONE_REFERENCE                       } from '../../../subworkflows/samtools/prepare.one.ref'
include { VCF_INPUT                                   } from '../../../subworkflows/nf/vcf_input'
include { DOWNLOAD_GTF_OR_GFF3  as  DOWNLOAD_GTF      } from '../../../modules/gtf/download'
include { DOWNLOAD_GNOMAD_SV                          } from '../../../modules/gnomad_sv/download.vcf'
include { GTF_TO_BED                                  } from '../../../modules/jvarkit/gtf2bed'
include { JVARKIT_VCFGNOMADSV                         } from '../../../modules/jvarkit/vcfgnomadsv'
include { BEDTOOLS_SLOP                               } from '../../../modules/bedtools/slop'
include { BEDTOOLS_MERGE                              } from '../../../modules/bedtools/merge'
include { BCFTOOLS_VIEW                               } from '../../../modules/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_INV          } from '../../../modules/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_INDEL        } from '../../../modules/bcftools/view'
include { ANNOTSV                                     } from '../../../subworkflows/annotsv'
include { DGV_DOWNLOAD                                } from '../../../modules/dgv/download'
include { READ_SAMPLESHEET                            } from '../../../subworkflows/nf/read_samplesheet'
include { META_TO_BAMS                                } from '../../../subworkflows/samtools/meta2bams1'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_INV        } from '../../../modules/bcftools/query'
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_INDEL      } from '../../../modules/bcftools/query'
include { BCFTOOLS_INDEX                              } from '../../../modules/bcftools/index'
include { COVERAGE_GRID                               } from '../../../modules/jvarkit/coveragegrid'
include { JVARKIT_SVLEN                               } from '../../../modules/jvarkit/svlen'
include { VCF_BED as JVARKIT_VCFBED_DGV               } from '../../../modules/jvarkit/vcfbed'
include { GHOSTSCRIPT_MERGE                           } from '../../../modules/gs/merge/main.nf'
include { PDF_NAVIGATION                              } from '../../../modules/pdf/navigation/main.nf'
include { GTF_ANNOTATION                              } from '../../../modules/gtf/annot1/main.nf'
include { PLOT_COVERAGE_01                            } from '../../../subworkflows/plotdepth'
include { BEDTOOLS_INTERSECT as INTERSECT_INV         } from '../../../modules/bedtools/intersect'
include { BEDTOOLS_INTERSECT as INTERSECT_INDEL         } from '../../../modules/bedtools/intersect'



workflow {
		validateParameters()

		if( params.help ) {
			log.info(paramsHelp())
			exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
		throw new IllegalArgumentException("undefined --samplesheet");
		}
	multiqc = Channel.empty()
	versions = Channel.empty()
	def metadata = [id:"annotsv"]
		

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

	META_TO_BAMS(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		READ_SAMPLESHEET.out.samplesheet
		)
	versions = versions.mix(META_TO_BAMS.out.versions)

	DOWNLOAD_GTF(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(DOWNLOAD_GTF.out.versions)

    DOWNLOAD_GNOMAD_SV( PREPARE_ONE_REFERENCE.out.dict )
    versions = versions.mix(DOWNLOAD_GNOMAD_SV.out.versions)

    DGV_DOWNLOAD( PREPARE_ONE_REFERENCE.out.dict )
    versions = versions.mix(DGV_DOWNLOAD.out.versions)


	VCF_INPUT(
        metadata.plus([
			path:params.vcf,
			require_index:true,
			arg_name:"vcf",
			required: true,
			unique : false
			])
        )
    versions = versions.mix(VCF_INPUT.out.versions)

	if(params.bed==null) {
		bed = PREPARE_ONE_REFERENCE.out.scatter_bed
	} else if(params.bed.matches("(coding[_-])?(exon|gene)s?")) {
		GTF_TO_BED(
			PREPARE_ONE_REFERENCE.out.dict,
			DOWNLOAD_GTF.out.gtf.map{meta,gtf,_tbi->[meta,gtf]}
			)
		versions = versions.mix(GTF_TO_BED.out.versions)
		bed =  GTF_TO_BED.out.bed
	} else {
		bed =Channel.of(params.bed).map{[[id:"bed"],file(f)]}
	}

	BEDTOOLS_SLOP(
		PREPARE_ONE_REFERENCE.out.fai,
		bed
		)
	versions = versions.mix(BEDTOOLS_SLOP.out.versions)

	BEDTOOLS_MERGE( BEDTOOLS_SLOP.out.bed )
	versions = versions.mix(BEDTOOLS_MERGE.out.versions)

	/** add SVLEN attribute if it is missing */
	JVARKIT_SVLEN(
		VCF_INPUT.out.vcf.map{meta,f,tbi->[meta,f]}
		)
	versions = versions.mix(JVARKIT_SVLEN.out.versions)

	/** index vcf */
	BCFTOOLS_INDEX(JVARKIT_SVLEN.out.vcf)
	versions = versions.mix(BCFTOOLS_INDEX.out.versions)

	/** remove small, large indels , non-PASS */
	BCFTOOLS_VIEW(
		BEDTOOLS_MERGE.out.bed,
		[[id:"no_sample"],[]],
		BCFTOOLS_INDEX.out.vcf
		)
	versions = versions.mix(BCFTOOLS_VIEW.out.versions)

	/** annotate with DGV */
	JVARKIT_VCFBED_DGV(
		DGV_DOWNLOAD.out.bed,
		BCFTOOLS_VIEW.out.vcf
		)
	versions = versions.mix(JVARKIT_VCFBED_DGV.out.versions)

	/** filter out frequent in gnomad */
	JVARKIT_VCFGNOMADSV(
		DOWNLOAD_GNOMAD_SV.out.vcf,
		JVARKIT_VCFBED_DGV.out.vcf
		)
	versions = versions.mix(JVARKIT_VCFGNOMADSV.out.versions)

	if(params.genes_of_interest==null) {
		genes_roi_ch = Channel.of([[id:"roi"],[]])
	} else {
		genes_roi_ch = Channel.of([[id:"roi"],file(params.genes_of_interest)])
	}

	GTF_ANNOTATION(
		PREPARE_ONE_REFERENCE.out.fai,
		DOWNLOAD_GTF.out.gtf.map{meta,gtf,tbi->[meta,gtf]},
		genes_roi_ch,
		JVARKIT_VCFGNOMADSV.out.vcf
		)
	versions = versions.mix(GTF_ANNOTATION.out.versions)


	ANNOTSV(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		GTF_ANNOTATION.out.vcf
		)
	versions = versions.mix(ANNOTSV.out.versions)

	BCFTOOLS_VIEW_INDEL(
		[[id:"nobed"],[]],
		[[id:"no_sample"],[]],
		GTF_ANNOTATION.out.vcf.map{meta,vcf->[meta,vcf,[]]}
		)
	versions = versions.mix(BCFTOOLS_VIEW_INDEL.out.versions)

	/** convert indels to BED */
	BCFTOOLS_QUERY_INDEL(BCFTOOLS_VIEW_INDEL.out.vcf)
	versions = versions.mix(BCFTOOLS_QUERY_INDEL.out.versions)
	
	if(params.genes_of_interest==null) {
		plot_indel_bed = BCFTOOLS_QUERY_INDEL.out.output
		}
	else
		{
		INTERSECT_INDEL(
			PREPARE_ONE_REFERENCE.out.fai,
				BCFTOOLS_QUERY_INDEL.out.output.combine(GTF_ANNOTATION.out.roi_bed)
					.map{meta1,a,meta2,b->[meta2,a,b]}
				)
		versions = versions.mix(INTERSECT_INDEL.out.versions)
		plot_indel_bed = INTERSECT_INDEL.out.bed
		}

	PLOT_COVERAGE_01(
		metadata,
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		DOWNLOAD_GTF.out.gtf,
		plot_indel_bed,
		META_TO_BAMS.out.bams
		)
	versions = versions.mix(PLOT_COVERAGE_01.out.versions)
	multiqc = multiqc.mix(PLOT_COVERAGE_01.out.multiqc)

	BCFTOOLS_VIEW_INV(
		[[id:"nobed"],[]],
		[[id:"no_sample"],[]],
		GTF_ANNOTATION.out.vcf.map{meta,vcf->[meta,vcf,[]]}
		)
	versions = versions.mix(BCFTOOLS_VIEW_INV.out.versions)



	MAKE_SAMPLESHEET2(
		META_TO_BAMS.out.bams
			.map{meta,bam,bai->"${meta.id}\tBAMS/${bam.name}\t${meta.status=="case"?"red":"green"}"}
			.collect()
			.map{L->[[id:"bams"],L.sort()]}
	)
	
	BCFTOOLS_QUERY_INV(BCFTOOLS_VIEW_INV.out.vcf)
	versions = versions.mix(BCFTOOLS_QUERY_INV.out.versions)
	
	INTERSECT_INV(
		PREPARE_ONE_REFERENCE.out.fai,
		BCFTOOLS_QUERY_INV.out.output.combine(GTF_ANNOTATION.out.roi_bed)
			.map{meta1,a,meta2,b->[meta2,a,b]}
		)
	versions = versions.mix(INTERSECT_INV.out.versions)

	COVERAGE_GRID(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		DOWNLOAD_GTF.out.gtf,
		DOWNLOAD_GNOMAD_SV.out.vcf,
		MAKE_SAMPLESHEET2.out.samplesheet,
		META_TO_BAMS.out.bams
			.map{meta,bam,bai->[bam,bai]}
			.flatMap()
			.collect()
			.map{L->[[id:"bams"],L.sort()]},
		BCFTOOLS_QUERY_INV.out.output
			.map{_meta,f->f}
			.splitCsv(header:false,sep:'\t')
			.map{row->[[id:"${row[0]}_${row[1]}_${row[2]}",title:"${row[0]}:${row[1]}-${row[2]}"],row[0],(row[1] as int),(row[2] as int)]}
		)
	versions = versions.mix(COVERAGE_GRID.out.versions)

	GHOSTSCRIPT_MERGE(COVERAGE_GRID.out.postscript)
	versions = versions.mix(GHOSTSCRIPT_MERGE.out.versions)

	PDF_NAVIGATION(GHOSTSCRIPT_MERGE.out.pdf.map{m,pdf->pdf}.collect().map{L->[[id:"inv"],L.sort()]})
	versions = versions.mix(PDF_NAVIGATION.out.versions)
}


process MAKE_SAMPLESHEET2 {
	executor "local"
	input:
		tuple val(meta),val(L)
	output:
		tuple val(meta),path("samplesheet.tsv"),emit:samplesheet
	script:
"""
cat << EOF > samplesheet.tsv
sample\tbam\tcolor
${L.join("\n")}
EOF
"""
}