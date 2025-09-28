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
nextflow.enable.dsl=2
include {assertKeyExistsAndNotEmpty          } from '../../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty            } from '../../../modules/utils/functions.nf'
include {assertKeyMatchRegex                 } from '../../../modules/utils/functions.nf'
include {VCF_STATS                           } from '../../../subworkflows/vcfstats'
include {COMPILE_VERSIONS                    } from '../../../modules/versions/main.nf'
include {MULTIQC                             } from '../../../modules/multiqc'
include {SCATTER_TO_BED                      } from '../../../subworkflows/gatk/scatterintervals2bed'
include {BEDTOOLS_MAKEWINDOWS                } from '../../../modules/bedtools/makewindows'
include {BED_CLUSTER                         } from '../../../modules/jvarkit/bedcluster'
include {GATK_BAM2VCF                        } from '../../../subworkflows/gatk/bam2vcf'
include {runOnComplete; dumpParams           } from '../../../modules/utils/functions.nf'
include {BCFTOOLS_GUESS_PLOIDY               } from '../../../modules/bcftools/guess_ploidy'
include {BCFTOOLS_STATS                      } from '../../../modules/bcftools/stats'


// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline ")
   exit 0
}


workflow {

	bams_ch = Channel.empty()
	versions = Channel.empty()
	multiqc = Channel.empty()

	ref_hash = [
		id:   file(params.fasta).simpleName,
		name: file(params.fasta).simpleName,
		ucsc_name : params.ucsc_name?:"undefined"
		]

	fasta = [ref_hash,file(params.fasta)]
	fai   = [ref_hash,file(params.fai)]
	dict  = [ref_hash,file(params.dict)]

    if(params.dbsnp!=null)
            {
            dbsnp =     [ref_hash, [] ,[] ]
            }
    else
            {
            dbsnp =     [ ref_hash, file(params.dbsnp), file(params.dbsnp+".tbi") ]
            }

    def gtf     = [ref_hash,[],[]]
  
    if(params.gtf!=null) {
        gtf = [ref_hash,file(params.gtf),file(params.gtf+".tbi")]
    }



	bams_and_ref = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')

    bams_ch = bams_and_ref
            .map{if(!it.containsKey("bam")) throw new IllegalArgumentException("${it} : bam missing"); return it;}
			.map{if(!(it.bam.endsWith(".cram") || it.bam.endsWith(".bam"))) throw new IllegalArgumentException("${it}.bam should end with bam or cram"); return it;}
			.map{it.bai?it: (it.bam.endsWith(".bam") ? it.plus(["bai":it.bam+".bai"]):  it.plus(["bai":it.bam+".crai"]))}
            .map{it.sample ? it : it.plus(["sample":it.bam.md5()])}
            .map{[[id:it.sample], file(it.bam),file(it.bai)]}
            .take(10)
            

    references = Channel.of([file(params.fasta),file(params.fai),file(params.dict)])

    references = references.mix(
            bams_and_ref
            .filter{it.containsKey("fasta")}
            .filter{testKeyExistsAndNotEmpty(it,"fasta")}
            .map{assertKeyMatchRegex(it,"fasta",".*\\.(fasta|fa|fna)")}
           
            .map{it.fai ? it : it.plus(["fai":it.fasta+".fai"])}
            .map{it.dict?it : it.plus(["dict":it.fasta.replaceAll("\\.(fasta|fa|fna)\$",".dict")])}
            .map{[file(it.fasta),file(it.fai),file(it.dict)]}
        )
    
    all_references = references
        .flatMap()
        .unique()
        .collect()
        .map{[ref_hash,it]}

    if(params.bed==null) {
        SCATTER_TO_BED(ref_hash,fasta,fai,dict)
        versions = versions.mix(SCATTER_TO_BED.out.versions)
        bed = SCATTER_TO_BED.out.bed
    } else {
        bed = Channel.of([ref_hash, file(params.bed)])
    }
  
   /* cut the bed/genome into parts for SV calling per region */
   BEDTOOLS_MAKEWINDOWS(bed)
   versions = versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

   /* if it's an exome , group the small genome together in BED */
   BED_CLUSTER(fasta,fai,dict,BEDTOOLS_MAKEWINDOWS.out.bed)
   versions = versions.mix(BED_CLUSTER.out.versions)
   beds_ch = BED_CLUSTER.out.bed
    .map{it[1]}
    .map{it instanceof List?it:[it]}
    .flatMap()
    .map{[[id:it.baseName],it]}
    .take(10)//TODO



    

GATK_BAM2VCF(
    ref_hash,
    fasta,
    fai,
    dict,
    dbsnp,
    all_references, //[meta, [ref files fa fai dict...]] all known reference
    beds_ch, // [meta,bed]
    bams_ch, // [meta,bam,bai]
    )
 versions = versions.mix(GATK_BAM2VCF.out.versions)
 
 /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
  BCFTOOLS_GUESS_PLOIDY(fasta, fai,GATK_BAM2VCF.out.vcf)
  versions = versions.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)



  /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
  BCFTOOLS_STATS(
    fasta,
    fai,
    bed,
    Channel.of(gtf).map{[it[0],it[1]]}.first(),//meta,gtf
    [[:],[]],//samples,
    GATK_BAM2VCF.out.vcf.map{[it[0],[it[1],it[2]]]}
    )
  versions = versions.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)


    COMPILE_VERSIONS(versions.collect().map{it.sort()})
    multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MULTIQC(multiqc.map{it[1]}.collect().map{[[id:"hapcaller"],it]})

}

runOnComplete(workflow)
