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
include { BCFTOOLS_CONCAT_CONTIGS                  } from '../../../subworkflows/bcftools/concat.contigs'
include { GTF_INPUT                                } from '../../../subworkflows/nf/gtf_input'
include { runOnComplete                            } from '../../../modules/utils/functions.nf'
include { parseBoolean                             } from '../../../modules/utils/functions.nf'
include { makeKey                                  } from '../../../modules/utils/functions.nf'
include { isBlank                                  } from '../../../modules/utils/functions.nf'
include { flatMapByIndex                           } from '../../../modules/utils/functions.nf'
include { PLINK_BFILE2VCF                          } from '../../../modules/plink/bfile2vcf'
include { REGENIE_STEP2                            } from '../../../modules/regenie/step2'
include { PREPARE_USER_BED                         } from '../../../subworkflows/bedtools/prepare.user.bed'
include { VCF_INPUT                                } from '../../../subworkflows/nf/vcf_input'
include { PREPARE_ONE_REFERENCE                    } from '../../../subworkflows/samtools/prepare.one.ref'
include { VCF_TO_CONTIGS                           } from '../../../subworkflows/bcftools/vcf2contigs'
include { DIGEST_SAMPLESHEET                       } from './sub.nf'
include { MDS_TO_COVARIATES                        } from './sub.nf'
include { UPDATE_PGEN                              } from './sub.nf'
include { FUNCTIONAL_ANNOTATION_SCORES             } from './sub.nf'
include { MAKE_FUNCTIONAL_ANNOT_PER_CTG            } from './sub.nf'
include { MERGE_REGENIE                            } from './sub.nf'
include { PLINK2_VCF2PGEN                          } from '../../../modules/plink/vcf2pgen'
include { PLINK2_MERGE_PGEN                        } from '../../../modules/plink/merge_pgen'
include { REGENIE_FUNCTIONAL_ANNOT                 } from '../../../modules/jvarkit/regeniefunctionalannot'
include { REGENIE_SLIDING_ANNOT                    } from '../../../modules/jvarkit/regenieslidingannot'
include { REGENIE_BED_ANNOT                        } from '../../../modules/jvarkit/regeniebedannot'
include { REGENIE_MAKE_ANNOT                       } from '../../../subworkflows/jvarkit/regeniemakeannot'
include { QQMAN                                    } from '../../../modules/qqman/main.nf'
include { ZIP                                      } from '../../../modules/utils/zip'
include { MULTIQC                                  } from '../../../modules/multiqc'
include { COMPILE_VERSIONS                         } from '../../../modules/versions'
include { JVARKIT_VCFSTATS                         } from '../../../modules/jvarkit/vcfstats'

workflow {
	versions = Channel.empty()
	multiqc = Channel.empty()
	metadata = [id:"regenie"]

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
	if(params.step1_loco==null) {
		log.warn("--step1_loco missing ")
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

	/** convert vcfs to PGEN */
	PLINK2_VCF2PGEN(
		PREPARE_ONE_REFERENCE.out.dict, //used for PAR regions
		DIGEST_SAMPLESHEET.out.plink_ped,
		VCF_TO_CONTIGS.out.vcf
			.filter{meta,_vcf,_tbi->parseBoolean(params.skip_XY)==false || !meta.contig.matches("(chr)?[XY]")}
			.map{meta,vcf,tbi->[meta, vcf]}
		)
	versions = versions.mix(PLINK2_VCF2PGEN.out.versions)

	/** merge all pgen files */
	PLINK2_MERGE_PGEN(
		PLINK2_VCF2PGEN.out.pgen
			.map{meta,pgen,psam,pvar->[pgen,psam,pvar]}
			.flatMap()
			.collect()
			.map{files->[metadata,files.sort()]}
		)
	versions = versions.mix(PLINK2_MERGE_PGEN.out.versions)

	keep_markers_ch = [[id:"pvar"],(params.step1_markers==null?([]):file("${params.step1_markers}"))]

	/** update pgen with data from pedigee */
	UPDATE_PGEN(
		DIGEST_SAMPLESHEET.out.plink_ped,
		PLINK2_MERGE_PGEN.out.pgen
		)
	versions = versions.mix(UPDATE_PGEN.out.versions)


	 /**
	 * COVARIATES, convert from MDS if needed
	 */
	if(params.covariates==null) {
		log.error("--covariates undefined")
		exit -1;
		}
	else if(params.covariates.endsWith(".mds")) {
		MDS_TO_COVARIATES([[id:"covariates"],file(params.covariates)])
		versions = versions.mix(MDS_TO_COVARIATES.out.versions)
		covar_ch = MDS_TO_COVARIATES.out.covariates
	} else {
		covar_ch = [[id:"covariates"],file(params.covariates)]
	}


		
	

	/**
	 * ANNOTATIONS
	 */
	the_annot_ch = Channel.empty()
	make_annot_ch = Channel.empty()

	if(parseBoolean(params.with_functional_annotations)) {
		
		if(params.custom_annotation2mask==null) {
			log.info("using default masks: ${moduleDir}/default_scores.txt")
			mask_ch = [[id:"masks"],file("${moduleDir}/default_scores.txt")]
			}
		else
			{
			mask_ch = [[id:"masks"],file("${params.custom_annotation2mask}")]
			}

		FUNCTIONAL_ANNOTATION_SCORES(mask_ch)
		versions = versions.mix(FUNCTIONAL_ANNOTATION_SCORES.out.versions)

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
			GTF_INPUT.out.gtf,
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
		make_annot_ch = make_annot_ch.mix(MAKE_SLIDING.out.tsv)
		versions = versions.mix(MAKE_SLIDING.out.versions)
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
		REGENIE_BED_ANNOT(
			ch1,
			VCF_TO_CONTIGS.out.vcf
				.filter{meta,_vcf,_tbi->parseBoolean(params.skip_XY)==false || !meta.contig.matches("(chr)?[XY]")}
				.map{meta,vcf,tbi->[meta.plus(id:meta.id+"."+meta.contig),vcf,tbi]}
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


	MERGE_REGENIE(
		PREPARE_ONE_REFERENCE.out.dict,
		REGENIE_STEP2.out.regenie
			.map{_meta,r->(r instanceof List?r:[r])}
			.flatMap()
			.collect()
			.map{files->[[id:"regenie"],files.sort()]}
		)
	versions = versions.mix(MERGE_REGENIE.out.versions)

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
