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

include { validateParameters                       } from 'plugin/nf-schema'
include { paramsHelp                               } from 'plugin/nf-schema'
include { paramsSummaryLog                         } from 'plugin/nf-schema'
include { samplesheetToList                        } from 'plugin/nf-schema'
include {ULTRA_RARES_ITERATION as ITER_10          } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_100         } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_250         } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_500         } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_1000        } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_ALL         } from './iteration.part.nf'
/*include {ULTRA_RARES_ITERATION as ITER_SECOND      } from './iteration.part.nf'
include {ULTRA_RARES_ITERATION as ITER_THIRD       } from './iteration.part.nf'*/
include { UNMAPPED                                 } from '../../subworkflows/unmapped'
include { runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { isBlank                                  } from '../../modules/utils/functions.nf'
include { SAMTOOLS_STATS                           } from '../../modules/samtools/stats'
include { MULTIQC                                  } from '../../modules/multiqc'
include { COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { META_TO_BAMS                             } from '../../subworkflows/samtools/meta2bams1'
include { READ_SAMPLESHEET                         } from '../../subworkflows/nf/read_samplesheet'
include { GTF_INPUT                                } from '../../subworkflows/nf/gtf_input'
include { GTF_TO_EXOME                             } from '../../modules/gtf/gtf2exome1'
include { BEDTOOLS_MAKEWINDOWS                     } from '../../modules/bedtools/makewindows'
include { BED_CLUSTER                              } from '../../modules/jvarkit/bedcluster'
include { SNPEFF_DOWNLOAD                          } from '../../modules/snpeff/download'
include { JVARKIT_VCF_TO_TABLE as VCF_TO_HTML      } from '../../modules/jvarkit/vcf2table'
include { JVARKIT_VCF_TO_TABLE as VCF_TO_TXT       } from '../../modules/jvarkit/vcf2table'
include { JVARKIT_FILTER_LOWQUAL                   } from '../../modules/jvarkit/lowqual'
include { GROUP_BY_GENE                            } from '../../modules/jvarkit/groupbygene'
include { META_TO_PED                              } from '../../subworkflows/pedigree/meta2ped'
include { VCF_INPUT as GNOMAD_INPUT                } from '../../subworkflows/nf/vcf_input'
include { ENCODE_BLACKLIST                         } from '../../modules/encode/blacklist' 
include { BEDTOOLS_SUBTRACT                        } from '../../modules/bedtools/subtract' 
include { BEDTOOLS_INTERSECT                       } from '../../modules/bedtools/intersect' 
include { BCFTOOLS_QUERY                           } from '../../modules/bcftools/query'
include { BCFTOOLS_SORT as BCFTOOLS_SORT2          } from '../../modules/bcftools/sort' 
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX2        } from '../../modules/bcftools/index' 
include { DOWNLOAD_CYTOBAND                        } from '../../modules/ucsc/download.cytobands'
include { DOWNLOAD_REFGENE                         } from '../../modules/ucsc/download.refgene'
include { IGV_REPORT                               } from '../../modules/igv/igv_report1'

List fractionate(List Lsrc,int count) {
	try {
		def L0 = Lsrc.collect() //clone
		def L = [];
		int i = 0;
		while(i< L0.size()) {
			def item = L0[i];
			if(!isBlank(item[0].status) && item[0].status=="case") {
				L.add(item);
				L0.remove(i);
				}
			else
				{
				++i;
				}
			}
		// sort using biggest files first, exomes at the end
		Collections.sort(L0,(A,B)->Long.compare(B[1].size(),A[1].size()))
		while((L.size()< count || count< 0) && !L0.isEmpty()) {
			L.add(L0.remove(0));
			}
		return L;
		}
	catch(Throwable err) {
		log.warn(String.valueOf(err.getMessage()));
		err.printStackTrace();
		throw err;
		}
	}


workflow {
	if(!workflow.stubRun) {
		validateParameters()
		}

	if( params.help ) {
		log.info(paramsHelp())
		exit 0
		}  else {
		// Print summary of supplied parameters
		log.info paramsSummaryLog(workflow)
		}



	def metadata = [id:"ultrares"]
	versions = Channel.empty()
	multiqc  = Channel.empty()


	if(params.fasta==null) {
			throw new IllegalArgumentException("undefined --fasta");
			}
	if(params.samplesheet==null) {
			throw new IllegalArgumentException("undefined --samplesheet");
			}

  /***************************************************
   *
   *  PREPARE FASTA REFERENCE
   *
   */
	PREPARE_ONE_REFERENCE(
			metadata,
			Channel.of(params.fasta).map{file(it)}.map{[[id:it.baseName],it]}
			)
	versions = versions.mix(PREPARE_ONE_REFERENCE.out.versions)


	/* no fastq samplesheet */
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

	/***************************************************
	*
	*  MAKE PEDIGREE
	*
	*/
	META_TO_PED(
		metadata,
		META_TO_BAMS.out.bams.map{meta,bam,bai->meta}
		)
	versions = versions.mix(META_TO_PED.out.versions)

	/***************************************************
	*
	*  MAKE BED
	*
	*/
	if(params.bed==null) {
		bed = PREPARE_ONE_REFERENCE.out.scatter_bed
		}
	else if(params.bed=="exome") {
		/***************************************************
		*
		*  DOWNLOAD GTF
		*
		*/
		GTF_INPUT(
				metadata.plus([
					arg_name: "gtf",
					require_index: true,
					download: true,
					path: params.gtf
					]),
				PREPARE_ONE_REFERENCE.out.dict
				)
		versions = versions.mix(GTF_INPUT.out.versions)

		GTF_TO_EXOME(
			PREPARE_ONE_REFERENCE.out.fai,
			GTF_INPUT.out.gtf.map{meta,gtf,tbi->[meta,gtf]}
			)
		versions = versions.mix(GTF_TO_EXOME.out.versions)

		bed  = GTF_TO_EXOME.out.bed
		}
	else
		{
		
		BEDTOOLS_INTERSECT(
			PREPARE_ONE_REFERENCE.out.fai,
			 Channel.of(params.bed)
			 	.map{bed->file(bed)}
				.map{bed->[[id:bed.baseName],bed]}
				.combine(PREPARE_ONE_REFERENCE.out.scatter_bed)
				.map{meta1,bed1,meta2,bed2->[meta1,bed1,bed2]}
			)
		versions = versions.mix(BEDTOOLS_INTERSECT.out.versions)
		bed = BEDTOOLS_INTERSECT.out.bed
		}

	/***************************************************
	*
	* Download encode black list
	*
	*/
	ENCODE_BLACKLIST(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(ENCODE_BLACKLIST.out.versions)
	BEDTOOLS_SUBTRACT(bed.combine(ENCODE_BLACKLIST.out.bed).map{meta1,bed1,_meta2,bed2->[meta1,bed1,bed2]})
    versions = versions.mix(BEDTOOLS_SUBTRACT.out.versions)
    bed = BEDTOOLS_SUBTRACT.out.bed.first()


	BEDTOOLS_MAKEWINDOWS( bed )
  	versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

	/* if it's an exome , group the small genome together in BED */
	BED_CLUSTER(
		PREPARE_ONE_REFERENCE.out.dict,
		BEDTOOLS_MAKEWINDOWS.out.bed
		)
	versions = versions.mix(BED_CLUSTER.out.versions)
	bed = BED_CLUSTER.out.bed
		.map{_meta,bed->bed}
		.map{it instanceof List?it:[it]}
		.flatMap()
		.map{bed->[[id:bed.baseName],bed]}
	

	bams_list10   = META_TO_BAMS.out.bams.collect(flat:false)
		.flatMap{fractionate(it,10)}
		.flatMap{meta,bam,bai->[bam,bai]}
		.collect(flat:true,sort:true)
		.map{files->[[id:"list10"],files]}

	bams_list100   = META_TO_BAMS.out.bams.collect(flat:false)
		.flatMap{fractionate(it,100)}
		.flatMap{meta,bam,bai->[bam,bai]}
		.collect(flat:true,sort:true)
		.map{files->[[id:"list100"],files]}

	bams_list250   = META_TO_BAMS.out.bams.collect(flat:false)
		.flatMap{fractionate(it,250)}
		.flatMap{meta,bam,bai->[bam,bai]}
		.collect(flat:true,sort:true)
		.map{files->[[id:"list250"],files]}

	bams_list500   = META_TO_BAMS.out.bams.collect(flat:false)
		.flatMap{fractionate(it,500)}
		.flatMap{meta,bam,bai->[bam,bai]}
		.collect(flat:true,sort:true)
		.map{files->[[id:"list500"],files]}
	
	bams_list1000   = META_TO_BAMS.out.bams.collect(flat:false)
		.flatMap{fractionate(it,1000)}
		.flatMap{meta,bam,bai->[bam,bai]}
		.collect(flat:true,sort:true)
		.map{files->[[id:"list1000"],files]}

	bams_listALL   = META_TO_BAMS.out.bams.collect(flat:false)
		.flatMap{fractionate(it,-1)}
		.flatMap{meta,bam,bai->[bam,bai]}
		.collect(flat:true,sort:true)
		.map{files->[[id:"list_all"],files]}

	/***************************************************
	*
	*  LOCATE GNOMAD
	*
	*/
    GNOMAD_INPUT(
      metadata.plus([
        arg_name : "gnomad",
        path: params.gnomad,
        require_index : true,
        required : true,
        unique : true
        ])
      )
    versions = versions.mix(GNOMAD_INPUT.out.versions)
	/***************************************************
	*
	*  FILTER JVARKIT
	*
	*/
	if(params.jvarkit_filter!=null) {
		jvarkit_filter = Channel.of(params.jvarkit_filter)
			.map{file(it)}
			.map{[[id:it.baseName],it]}
			.first()
		}
	else
		{
		jvarkit_filter = META_TO_PED.out.select_code
		}
	
	/***************************************************
	*
	* DOWNLOAD SNPEFF
	*
	*/
	SNPEFF_DOWNLOAD(PREPARE_ONE_REFERENCE.out.fai
			.map{meta,_fai->[meta,file(params.snpeff_database_directory)]}
		)
	versions = versions.mix(SNPEFF_DOWNLOAD.out.versions)

	ITER_10(
		metadata.plus([id:"list10"]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		SNPEFF_DOWNLOAD.out.database,
		META_TO_PED.out.pedigree_gatk,
		GNOMAD_INPUT.out.vcf.first(),
		jvarkit_filter,
		bed,
		bams_list10
		)
	versions = versions.mix(ITER_10.out.versions)
	multiqc = multiqc.mix(ITER_10.out.multiqc)



	ITER_100(
		metadata.plus([id:"list100"]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		SNPEFF_DOWNLOAD.out.database,
		META_TO_PED.out.pedigree_gatk,
		GNOMAD_INPUT.out.vcf.first(),
		jvarkit_filter,
		ITER_10.out.bed,
		bams_list100
		)
	versions = versions.mix(ITER_100.out.versions)
	multiqc = multiqc.mix(ITER_100.out.multiqc)

	ITER_250(
		metadata.plus([id:"list250"]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		SNPEFF_DOWNLOAD.out.database,
		META_TO_PED.out.pedigree_gatk,
		GNOMAD_INPUT.out.vcf.first(),
		jvarkit_filter,
		ITER_100.out.bed,
		bams_list250
		)
	versions = versions.mix(ITER_250.out.versions)
	multiqc = multiqc.mix(ITER_250.out.multiqc)


	ITER_500(
		metadata.plus([id:"list500"]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		SNPEFF_DOWNLOAD.out.database,
		META_TO_PED.out.pedigree_gatk,
		GNOMAD_INPUT.out.vcf.first(),
		jvarkit_filter,
		ITER_250.out.bed,
		bams_list500
		)
	versions = versions.mix(ITER_500.out.versions)
	multiqc = multiqc.mix(ITER_500.out.multiqc)

	ITER_1000(
		metadata.plus([id:"list1000"]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		SNPEFF_DOWNLOAD.out.database,
		META_TO_PED.out.pedigree_gatk,
		GNOMAD_INPUT.out.vcf.first(),
		jvarkit_filter,
		ITER_500.out.bed,
		bams_list1000
		)
	versions = versions.mix(ITER_1000.out.versions)
	multiqc = multiqc.mix(ITER_1000.out.multiqc)

	ITER_ALL(
		metadata.plus([id:"list_all"]),
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		SNPEFF_DOWNLOAD.out.database,
		META_TO_PED.out.pedigree_gatk,
		GNOMAD_INPUT.out.vcf.first(),
		jvarkit_filter,
		ITER_1000.out.bed,
		bams_listALL
		)
	versions = versions.mix(ITER_ALL.out.versions)
	multiqc = multiqc.mix(ITER_ALL.out.multiqc)

	JVARKIT_FILTER_LOWQUAL(ITER_ALL.out.vcf.map{meta,vcf,tbi->[meta,vcf]})
	versions = versions.mix(JVARKIT_FILTER_LOWQUAL.out.versions)

	VCF_TO_HTML(
		META_TO_PED.out.pedigree,
		JVARKIT_FILTER_LOWQUAL.out.vcf
		)
	versions = versions.mix(VCF_TO_HTML.out.versions)
	VCF_TO_TXT(
		META_TO_PED.out.pedigree,
		JVARKIT_FILTER_LOWQUAL.out.vcf
		)
	versions = versions.mix(VCF_TO_TXT.out.versions)


	DOWNLOAD_CYTOBAND(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(DOWNLOAD_CYTOBAND.out.versions)
	DOWNLOAD_REFGENE(PREPARE_ONE_REFERENCE.out.dict)
	versions = versions.mix(DOWNLOAD_REFGENE.out.versions)

	BCFTOOLS_SORT2(JVARKIT_FILTER_LOWQUAL.out.vcf)
	versions = versions.mix(BCFTOOLS_SORT2.out.versions)

	BCFTOOLS_INDEX2(BCFTOOLS_SORT2.out.vcf)
	versions = versions.mix(BCFTOOLS_INDEX2.out.versions)


	GROUP_BY_GENE(
		META_TO_PED.out.cases,
		META_TO_PED.out.controls,
		BCFTOOLS_INDEX2.out.vcf
		)
	versions = versions.mix(GROUP_BY_GENE.out.versions)

	BCFTOOLS_QUERY(BCFTOOLS_INDEX2.out.vcf.map{meta,vcf,_tbi->[meta,vcf]})
	versions = versions.mix(BCFTOOLS_QUERY.out.versions)
	ch1 = BCFTOOLS_QUERY.out.output
		.map{meta,f->f}
		.splitCsv(header:false,sep:'\t')
		.map{v->[
			id : String.join("_",v[0],v[1],v[2],v[3]),
			contig: v[0],
			start : v[1],
			end : v[2],
			ref: v[3]
			]}
		.combine(META_TO_BAMS.out.bams
			.filter{meta,bam,bai->meta.status=="case"}
			.map{meta,bam,bai->[bam,bai]}
			)
		.groupTuple()


	

	IGV_REPORT(
		PREPARE_ONE_REFERENCE.out.fasta,
		PREPARE_ONE_REFERENCE.out.fai,
		PREPARE_ONE_REFERENCE.out.dict,
		DOWNLOAD_CYTOBAND.out.bed,
		DOWNLOAD_REFGENE.out.tabix,
		[[id:"no_sample"],[]],
		META_TO_PED.out.pedigree,
		BCFTOOLS_INDEX2.out.vcf,
		ch1
		)
	versions = versions.mix(IGV_REPORT.out.versions)


	README(
		metadata,
		jvarkit_filter
		)
}


process README {
executor "local"
input:
	val(meta)
	tuple val(meta2),path(select_code)
output:
	tuple val(meta),path("README.md"),emit:readme
script:
"""


if ${select_code?true:false}
then
cat << '__EOF__' >> README.md
## Variant Filtration

VCF was filtered with jvarkit:vcffilterjdk :
<pre>
__EOF__

cat "${select_code}" >> README.md

cat << '__EOF__' >> README.md
</pre>
__EOF__

fi



cat << '__EOF__' >> README.md

# Output

  - <code>*.list_all.genes.bed</code> contient la list des genes trouves.
  - <code>*.list_all.sort.vcf.gz</code> les variants trouves, sous forme de fichier VCF
  - <code>*.list_all.table.html</code> les variants trouvés, sous forme de table html à ouvrir dans un navigateur web
  - <code>*.list_all.table.txt</code> les variants trouvés, sous forme de table en texte.
  - <code>pipeline_info</code> meta data about the pipeline
  - <code>IGV/*</code> context of the variants viewed as web-IGV

<pre>
|-- ${params.prefix}.list_all.genes.bed
|-- ${params.prefix}.list_all.sort.vcf.gz
|-- ${params.prefix}.list_all.table.html
|-- ${params.prefix}.list_all.table.txt
|-- IGV
|   |-- chr1_9377032_9377032_A.igv.html
|   |-- chr1_130007156_130007156_C.igv.html
|   |-- chr1_20570747_20570747_G.igv.html
|   |-- chr1_56584558_56584558_G.igv.html
|   +-- chr1_36235043_36235043_G.igv.html
+-- pipeline_info
    +-- execution_report_2025-11-29_22-00-43.html
</pre>
__EOF__

"""
}
