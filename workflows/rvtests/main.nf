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
	 *
	 * NORMALIZE_VCF
	 *
	 */
	NORMALIZE_VCF(VCF_INPUT.out.vcf)
	versions = versions.mix(NORMALIZE_VCF.out.versions)
	vcf_ch = NORMALIZE_VCF.out.vcf
		.map{meta,vcf,tbi->[meta.plus(id:makeKey(vcf,meta.contig)),vcf,tbi]} // git a unique id vcf/contig to vcf stream


	set_file_ch = Channel.empty()

	/**
	 * ANNOTATIONS
	 */

	if(parseBoolean(params.with_functional_annotations)) {
		SETFILE_FOR_FUNCTIONAL_ANNOT(vcf_ch)
		versions = versions.mix(SETFILE_FOR_FUNCTIONAL_ANNOT.out.versions)
		
		set_file_ch = set_file_ch.mix(
			SETFILE_FOR_FUNCTIONAL_ANNOT.out.setfile
				.flatMap{row->flatMapByIndex(row,1)}
				.join(vcf_ch)
				.map{meta,setfile,vcf,tbi->[meta.plus(id:makeKey(setfile),type:"functional"),setfile,vcf,tbi]}
			)
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

		SETFILE_FOR_SLIDING_WINDOWS(
			vcf_ch
				.combine(windows_ch)
				.map{meta,vcf,tbi,w_size,w_shift->[meta.plus(win_size:w_size,win_shift:w_shift,id:meta.id+"_"+meta.contig),vcf,tbi]}
			)
		versions = versions.mix(SETFILE_FOR_SLIDING_WINDOWS.out.versions)

		set_file_ch = set_file_ch.mix()
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

		dispatch_ch = 
			.multiMap{meta1,bed,meta2,vcf,tbi->
				bed: [meta1,bed]
				vcf: [meta2.plus(id:makeKey(meta2.id+"."+meta2.contig+"."+meta1.id)),vcf,tbi]
				}


		SETFILE_FOR_BED(
			ch1.combine(
				vcf_ch
					.flatMap{meta,vcf,idx->[vcf,idx]}.
					.collect()
					.map{f->f.sort()}
				)
			)
		set_file_ch = set_file_ch.mix(SETFILE_FOR_BED.out.tsv)
		versions = versions.mix(SETFILE_FOR_BED.out.versions)
		}
	else
		{
		log.info("NO custom contig/start/end/annot/gene/file.");
		}

	/** 
	 * CONVERT TSV data to regenie input
	 *
	 */
	LINUX_SPLIT(set_file_ch)
	versions = versions.mix(LINUX_SPLIT.out.versions)

	
	RVTESTS_APPLY(
		pedigree,
		vcf,
		setFile
		)
	versions = versions.mix(RVTESTS_APPLY.out.versions)

	
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


	}

runOnComplete(workflow)


/**
 * normalize vcf for rvtests, changing chromosome notation 'chr1' to '1'
 */
process NORMALIZE_VCF {
label "process_single"
tag "${meta.id}"
input:
	tuple val(meta),path(vcf),path(tbi)
output:
	tuple val(meta),path(vcf),path(tbi),emit:vcf
script:
	def contig = meta.contig
	def prefix  = task.ext.prefix?:"${meta.id}.${contig}"
"""
echo "${contig}\t${contig}" | sed 's/\tchr/\t' > TMP/contigs.txt

bcftools annotate \\
	--threads ${task.cpus} \\
	--regions "${contig}"  \\
	${contig.startsWith("chr")?"--rename-chrs TMP/contigs.txt":""} \\
	-x 'FILTER,ID,QUAL,^FORMAT/GT' \\
	-O z \\
	-o TMP/jeter.vcf.gz

bcftools index --threads ${task.cpus} -f -t  TMP/jeter.vcf.gz

mv TMP/jeter.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter.vcf.gz.tbi ${prefix}.vcf.gz.tbi

touch versions.yml
"""
stub:
	def prefix  = task.ext.prefix?:"${meta.id}.${contig}"
"""
touch versions.yml ${prefix}.norm.vcf.gz ${prefix}.norm.vcf.gz.tbi
"""
}

process SETFILE_FOR_FUNCTIONAL_ANNOT {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(vcf),path(tbi)
output:
	tuple val(meta),path("*.setFile"),optional:true,emit:setFile
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}.func"
	def args1 = task.ext.args1?:""
"""
mkdir -p TMP/OUT1
bcftools view -G  ${args1} ${vcf} |\\
${jvarkit} vcffilterso --acn '' -r -R |\\
${jvarkit} vcfgenesplitter --extractors '' -o TMP/OUT1

find TMP/OUT1 -type f -name "*.vcf.gz" | sort -T TMP | uniq > TMP/jeter.list

cat TMP/jeter.list | while read F
do
	basename "\${F}" .vcf.gz | tr "\n" , >> TMP/jeter.setFile
	bcftools query -f '%CHROM:%POS-%END' |  paste -s -d',' >> TMP/jeter.setFile
done

mv TMP/jeter.setFile ${prefix}.setFile || true

touch versions.yml
"""
stub:
"""
touch versions.yml
"""
}

process SETFILE_FOR_SLIDING_WINDOWS {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(vcf),path(tbi)
output:
	tuple val(meta),path("*.setFile"),optional:true,emit:setFile
	path("versions.yml"),emit:versions
script:
	def contig = meta.contig
	def args1 = task.ext.args1?:""
	def win_size = meta.win_size
	def win_shift = meta.win_shift
	def prefix = task.ext.prefix?:"${meta.id}.sliding.w${win_size}_s${win_shift}"
"""
mkdir -p TMP/OUT1


# create sliding bed
bcftools index -s "${vcf}" |\\
	awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
	bedtools makewindows -w ${win_size} ${win_shift} |\\
	awk -F '\t' '{printf("%s:%d-%s\\n",\$1,int(\$2)+1,\$3);}' | while read RGN
	do
		bcftools query --regions "\${RGN}" ${args1} -f '%CHROM:%POS-%END'  ${vcf} |\\
			uniq |\\
			paste -s -d','|\\
			awk -v RGN="\${RGN}" 'BEGIN {R2=RGN;gsub(/[^A-Za-z0-9]+/,"_",R2);} {printf("%s.sliding.w${win_size}_s${win_shift}\t%s\\n",R2,\$0);}' >> TMP/jeter.setFile
	done

mv TMP/jeter.setFile ${prefix}.setFile || true

touch versions.yml
"""
stub:
"""
touch versions.yml 
"""
}


process SETFILE_FOR_BED {
label "process_single"
tag "${meta.id}"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(bed)
	tuple val(meta),path("VCF/*")
output:
	tuple val(meta),path("*.setFile"),optional:true,emit:setFile
	path("versions.yml"),emit:versions
script:
	def contig = meta.contig
	def args1 = task.ext.args1?:""
	def win_size = meta.win_size
	def win_shift = meta.win_shift
	def prefix = task.ext.prefix?:"${meta.id}.bed${meta1.id}"
"""
mkdir -p TMP/OUT1

find VCF -name "*.vcf.gz" > TMP/jeter.list



${bed.name.endsWith(".gz")?"gunzip -c":"cat"} "${bed}" |\\
	sed 's/^chr//' |\\
	awk -F '\t' '{S=\$4;if(NF<4 || S=="") {S=sprintf("%s_%s_%s",\$1,\$2,\$3);} printf("%s\t%s\t%s\\n",\$1,\$2,\$3,S);}' |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq > TMP/jeter1.bed

cut -f1,2,3  TMP/jeter1.bed |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	bedtools merge > TMP/jeter2.bed

bcftools concat --drop-genotypes ${args1} --file-list TMP/jeter.list --regions-file TMP/jeter2.bed  -O u |\\
	bcftools query --regions-file TMP/jeter.bed -f '%CHROM\t%POS0\t%END' |\\
	sort -T TMP -t '\t' -k1,1 -k2,2n |\\
	uniq > TMP/jeter3.bed

bedtools intersect -a TMP/jeter1.bed -b TMP/jeter2.bed -wa -wb |\\
	awk -F '\t' '{printf("%s\t%s:%s-%s\\n",\$4,\$1,int(\$2),\$3);}' |
	datamash groupby 1 unique 2 | while read RGNID RGN
	do
		echo "\${RGN}" | tr ":-" "\t" > TMP/jeter.bed
		bcf
		
			uniq |\\
			paste -s -d','|\\
			awk -v RGN="\${RGNID}" 'BEGIN {R2=RGNID;gsub(/[^A-Za-z0-9]+/,"_",R2);} {printf("\t%s\\n",R2,\$0);}' >> TMP/jeter.setFile
	done

mv TMP/jeter.setFile ${prefix}.setFile || true

touch versions.yml
"""
stub:
"""
touch versions.yml 
"""
}