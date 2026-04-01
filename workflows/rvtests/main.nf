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
/** https://github.com/FINNGEN/regenie-pipelines */
// https://gitlab.univ-nantes.fr/pierre.lindenbaum/gazoduc-nf/-/blob/d9f89be7c2cccc39336ea282bd54e26146ca645e/workflows/regenie/main.nf

include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include { BCFTOOLS_CONCAT_CONTIGS                  } from '../../subworkflows/bcftools/concat.contigs'
include { GTF_INPUT                                } from '../../subworkflows/nf/gtf_input'
include { runOnComplete                            } from '../../modules/utils/functions.nf'
include { parseBoolean                             } from '../../modules/utils/functions.nf'
include { makeKey                                  } from '../../modules/utils/functions.nf'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { flatMapByIndex                           } from '../../modules/utils/functions.nf'
include { PLINK_BFILE2VCF                          } from '../../modules/plink/bfile2vcf'
include { REGENIE_STEP2                            } from '../../modules/regenie/step2'
include { PREPARE_USER_BED                         } from '../../subworkflows/bedtools/prepare.user.bed'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { VCF_TO_CONTIGS                           } from '../../subworkflows/bcftools/vcf2contigs'

include { QQMAN                                    } from '../../modules/qqman/main.nf'
include { ZIP                                      } from '../../modules/utils/zip'
include { MULTIQC                                  } from '../../modules/multiqc'
include { COMPILE_VERSIONS                         } from '../../modules/versions'
include { JVARKIT_VCFSTATS                         } from '../../modules/jvarkit/vcfstats'
include { CIRCULAR_MANHATTAN                       } from '../../subworkflows/circular/manhattan'
include { BATIK                                    } from '../../subworkflows/batik'


workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()
	metadata = [id:"rvtests"]

	if(params.vcf==null) {
		log.warn("--vcf missing")
		exit -1
		}
	if(params.fasta==null) {
		log.warn("--fasta missing")
		exit -1
		}
	if(params.covariates==null) {
		log.warn("--covariates missing")
		exit -1
	}

	if(params.samplesheet==null) {
		log.warn("--samplesheet missing")
		exit -1
	}

  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata.plus(with_scatter:false),
			Channel.of(params.fasta).map{f->file(f)}.map{f->[[id:f.baseName],f]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


   /***************************************************
   *
   * PREPARE BED
   *
   */
  if(params.bed==null) {
		bed =  PREPARE_ONE_REFERENCE.out.scatter_bed
		}
	else 
		{
		PREPARE_USER_BED(
				metadata.plus([:]),
				PREPARE_ONE_REFERENCE.out.fasta,
				PREPARE_ONE_REFERENCE.out.fai,
				PREPARE_ONE_REFERENCE.out.dict,
				PREPARE_ONE_REFERENCE.out.scatter_bed,
				Channel.of([[id:"capture"],file(params.bed)])
				)
		versions = versions.mix(PREPARE_USER_BED.out.versions)
		multiqc = multiqc.mix(PREPARE_USER_BED.out.multiqc)
		bed = PREPARE_USER_BED.out.bed
		}

   /***************************************************
   *
   * VCF_INPUT
   *
   */
	VCF_INPUT(metadata.plus([
			path: params.vcf,
			arg_name: "vcf",
			require_index : true,
			required: true,
			unique : false
			]))
	versions = versions.mix(VCF_INPUT.out.versions)
	
	/** extract pedigrees, cases from vcf input and samplesheet */
	DIGEST_SAMPLESHEET (
		VCF_INPUT.out.vcf.take(1).map{meta,vcf,_tbi->[meta,vcf]},
		[[id:"samplesheet"],file(params.samplesheet)]
		)
	versions = versions.mix(DIGEST_SAMPLESHEET.out.versions)
	
	/** extract contig for each vcf */
	VCF_TO_CONTIGS(metadata, VCF_INPUT.out.vcf)
	versions = versions.mix(VCF_TO_CONTIGS.out.versions)

	/* check the only one vcf, or one contig/vcf */
	VCF_TO_CONTIGS.out.vcf
		.map{meta,vcf,_tbi->[meta.contig,vcf]}
		.groupTuple()
		.filter{contig,array->array.size()!=1}
		.map{contig,array->"ERROR  : multiple vcf for contig ${contig} e.g: ${array[0]} and ${array[1]}."}
		.mix(VCF_INPUT.out.vcf
			.count()
			.map{c->(c==1L?"":"More than one vcf")}
			.filter{s->!isBlank(s)}
			)
			.collect()
			.filter{array->array.size() > 1} /* two error flags */
			.subscribe{array->throw new IllegalArgumentException("${array.join(" ")}")}
	
	/*run statistics over the vcf */
	JVARKIT_VCFSTATS(
		DIGEST_SAMPLESHEET.out.sample2population,
		VCF_INPUT.out.vcf
			.map{meta,vcf,tbi->[[id:"input_vcf"],vcf,tbi]}
			.groupTuple()
			.map{meta,vcf,tbi->[meta,vcf.sort(),tbi.sort()]}
		)
	versions = versions.mix(JVARKIT_VCFSTATS.out.versions)
	multiqc = multiqc.mix(JVARKIT_VCFSTATS.out.json)


	

	/**
	 * ANNOTATIONS
	 */
	the_annot_ch = Channel.empty()
	make_annot_ch = Channel.empty()

	if(parseBoolean(params.with_functional_annotations)) {
		
		GTF_INPUT(
			metadata.plus(
				arg_name:"gtf",
				require_index:false,
				download:true,
				path: params.gtf
			),
			PREPARE_ONE_REFERENCE.out.dict
			)
		versions = versions.mix(GTF_INPUT.out.versions)

		REGENIE_FUNCTIONAL_ANNOT(
			FUNCTIONAL_ANNOTATION_SCORES.out.tsv,
			GTF_INPUT.out.gtf.map{meta,gtf,_tbi->[meta,gtf]},
			VCF_TO_CONTIGS.out.vcf
				.filter{meta,_vcf,_tbi->parseBoolean(params.skip_XY)==false || !meta.contig.matches("(chr)?[XY]")}
				.map{meta,vcf,tbi->[meta.plus(id:meta.id+".annot."+meta.contig),vcf,tbi]}
			)
		versions = versions.mix(REGENIE_FUNCTIONAL_ANNOT.out.versions)
		make_annot_ch = make_annot_ch.mix(REGENIE_FUNCTIONAL_ANNOT.out.tsv)

		//the_annot_ch = the_annot_ch.mix(REGENIE_MAKE_ANNOT.out.tsv)
		}

	/**
	 *
	 * Sliding windows for Regenie
	 * A pair of integers window_size,window_shift
	 *
	 */ 
	if(!isBlank(params.sliding_windows)) {
		windows_ch = Channel.of("${params.sliding_windows}").
			flatMap{T->{
			def L=[];
			def a = T.split("[,]");
			if(a.size()%2!=0) throw new IllegalArgumentException("not a pair ${params.sliding}");
			for(int i=0;i<a.size();i+=2) {
				L.add([ (a[i+0] as int) , (a[i+1] as int) ]);
				}
			return L;
			}}

		REGENIE_SLIDING_ANNOT(
			VCF_TO_CONTIGS.out.vcf
				.filter{meta,_vcf,_tbi->parseBoolean(params.skip_XY)==false || !meta.contig.matches("(chr)?[XY]")}
				.combine(windows_ch)
				.map{meta,vcf,tbi,w_size,w_shift->[meta.plus(win_size:w_size,win_shift:w_shift,id:meta.id+"_"+meta.contig),vcf,tbi]}
			)
		make_annot_ch = make_annot_ch.mix(REGENIE_SLIDING_ANNOT.out.tsv)
		versions = versions.mix(REGENIE_SLIDING_ANNOT.out.versions)
		}

	/**
	 *
	 * Custom BED file of annotations
	 *
	 */ 
	if(params.select_bed!=null) {
		log.info("making bed from ${params.select_bed}");
		if(params.select_bed.endsWith(".list")) {
			ch1 = Channel.fromPath(params.select_bed)
				.splitText()
				.map{fn->[[id:makeKey(fn)],file(nf.trim())]}
				;
			}
		else
			{
			ch1  = Channel.of(params.select_bed)
				.map{fn->[[id:makeKey(fn)],file(fn)]}
				;
			}

		dispatch_ch = ch1.combine(
				VCF_TO_CONTIGS.out.vcf
					.filter{meta,_vcf,_tbi->parseBoolean(params.skip_XY)==false || !meta.contig.matches("(chr)?[XY]")}
				)
			.multiMap{meta1,bed,meta2,vcf,tbi->
				bed: [meta1,bed]
				vcf: [meta2.plus(id:makeKey(meta2.id+"."+meta2.contig+"."+meta1.id)),vcf,tbi]
				}


		REGENIE_BED_ANNOT(
			dispatch_ch.bed,
			dispatch_ch.vcf
			)
		make_annot_ch = make_annot_ch.mix(REGENIE_BED_ANNOT.out.tsv)
		versions = versions.mix(REGENIE_BED_ANNOT.out.versions)
		}
	else
		{
		log.info("NO custom contig/start/end/annot/gene/file.");
		}

	/** 
	 * CONVERT TSV data to regenie input
	 *
	 */
	REGENIE_MAKE_ANNOT(
			metadata,
			FUNCTIONAL_ANNOTATION_SCORES.out.tsv,
			make_annot_ch
			)
	versions = versions.mix(REGENIE_MAKE_ANNOT.out.versions)
	multiqc = multiqc.mix(REGENIE_MAKE_ANNOT.out.multiqc)


	
	REGENIE_STEP2(
		UPDATE_PGEN.out.pgen.first(),
		covar_ch,
		DIGEST_SAMPLESHEET.out.plink_ped.first(), 
		[[id:"locostep1"],file(params.step1_loco)] ,
		REGENIE_MAKE_ANNOT.out.annotations
		)
	versions = versions.mix(REGENIE_STEP2.out.versions)

	/** merge all results and group all data by test/frequency */
	MERGE_REGENIE(
		PREPARE_ONE_REFERENCE.out.dict,
		REGENIE_STEP2.out.regenie
			.map{_meta,r->(r instanceof List?r:[r])}
			.flatMap()
			.collect()
			.map{files->[[id:"regenie"],files.sort()]}
		)
	versions = versions.mix(MERGE_REGENIE.out.versions)

	/** for multuqc, generate a table with the best hits */
	BEST_HITS(MERGE_REGENIE.out.tsv)
	versions = versions.mix(BEST_HITS.out.versions)
	multiqc = multiqc.mix(BEST_HITS.out.tsv)

	//join manifest data (containing test-name, freq) and the tsv for the same test
	to_qqman = MERGE_REGENIE.out.manifest
		.map{meta,mf->mf}
		.splitCsv(header:true,sep:'\t')
		.combine(
			MERGE_REGENIE.out.regenie
				.map{_meta,r->(r instanceof List?r:[r])}
				.flatMap()
			)
		.filter{meta,f->meta.filename==f.name}
		.map{meta,f->[meta.plus(id:meta.filename),f]}
	
	QQMAN(to_qqman)
	versions = versions.mix(QQMAN.out.versions)
	multiqc = multiqc.mix(QQMAN.out.manhattan)
	multiqc = multiqc.mix(QQMAN.out.qqplot)

	/**
	 * Make circos plots
	 */
	if(parseBoolean(params.with_circos)) {
		/** pdf can be big in size, just plot if LOG10P > 'x' */
		to_manhattan = to_qqman
			.map{meta,f->[meta,f,f]}/* duplicate regenie file */
			.splitCsv(header:true,sep:' '/*space */)
			.filter{_meta,row,_f->row.LOG10P!=null && row.LOG10P!="NA" && (row.LOG10P as double)>= (params.circos_treshold as double)} /* keep files having good p_value */
			.map{meta,_row,f->[meta,f]}
			.unique()


		CIRCULAR_MANHATTAN(
			metadata,
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			GTF_INPUT.out.gtf.map{meta,gtf,_tbi->[meta,gtf]},
			to_manhattan
			)
			
		versions = versions.mix(CIRCULAR_MANHATTAN.out.versions)
		multiqc = multiqc.mix(CIRCULAR_MANHATTAN.out.multiqc)

		/** convert SVG to pdf/png... */
		BATIK(
			metadata.plus(
				with_pdf:true,
				with_png:true,
				with_jpg:false
				),
			CIRCULAR_MANHATTAN.out.svg
			)
		versions = versions.mix(BATIK.out.versions)
		multiqc = multiqc.mix(BATIK.out.multiqc)
		}


	/**
	 * At the end, make QC, make versions, zip results
	 */
	COMPILE_VERSIONS(versions.collect().map{it.sort()})

	MULTIQC(
		[[id:"noconfig"],[]],
		multiqc.map{meta,f->f}
			.mix(COMPILE_VERSIONS.out.multiqc)
			.collect()
			.map{files->[[id:"regenie_mqc"],files.sort()]}
		)

	ZIP(
		QQMAN.out.manhattan
			.mix(QQMAN.out.qqplot)
			.mix(MERGE_REGENIE.out.tsv)
			.map{_meta,f->f}
			.collect()
			.map{files->[[id:"regenie"],files.sort()]}
		)
	versions = versions.mix(ZIP.out.versions)
	/*
	MERGE_AND_PLOT(
		reference,
		step2_ch.output.groupTuple()
		)
	
	ANNOT_HITS(reference,ch5.output.collect())


	
	MAKE_PNG_ARCHIVE(ch5.images.mix(ch5.ascii).flatten().collect())

	MAKE_SNPLIST(
		step2_ch.output.map{it[1]}.collect(),
		step2_ch.masks_snplist.map{it[1]}.collect()
		)

	README()
	*/

	}

runOnComplete(workflow)
