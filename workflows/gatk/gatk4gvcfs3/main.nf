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
include {assertKeyExistsAndNotEmpty          } from '../../../modules/utils/functions.nf'
include {testKeyExistsAndNotEmpty            } from '../../../modules/utils/functions.nf'
include {assertKeyMatchRegex                 } from '../../../modules/utils/functions.nf'
include {VCF_STATS                           } from '../../../subworkflows/vcfstats'
include {SCATTER_TO_BED                      } from '../../../subworkflows/gatk/scatterintervals2bed'
include {BEDTOOLS_MAKEWINDOWS                } from '../../../modules/bedtools/makewindows'
include {BED_CLUSTER                         } from '../../../modules/jvarkit/bedcluster'
include {GATK_BAM2VCF                        } from '../../../subworkflows/gatk/bam2vcf'
include {runOnComplete; dumpParams           } from '../../../modules/utils/functions.nf'
include {BCFTOOLS_GUESS_PLOIDY               } from '../../../modules/bcftools/guess_ploidy'
include {BCFTOOLS_STATS                      } from '../../../modules/bcftools/stats'
include {HAPLOTYPECALLER                     } from '../../../subworkflows/gatk/haplotypecaller'
include {HAPLOTYPECALLER_DIRECT              } from '../../../subworkflows/gatk/haplotypecaller.direct'
include {SAMTOOLS_SAMPLES                    } from '../../../modules/samtools/samples'
include {META_TO_PED                         } from '../../../subworkflows/pedigree/meta2ped'
include {PREPARE_REFERENCE                   } from '../../../subworkflows/samtools/prepare.ref'
include {MULTIQC                             } from '../../../subworkflows/multiqc'


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

    PREPARE_REFERENCE(ref_hash, fasta)
    versions = versions.mix(PREPARE_REFERENCE.out.versions)
	fai   = PREPARE_REFERENCE.out.fai
	dict  = PREPARE_REFERENCE.out.dict



    if(params.dbsnp==null)
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



  
    
    all_references =  bams_and_ref
        .filter{it.containsKey("fasta")}
        .filter{testKeyExistsAndNotEmpty(it,"fasta")}
        .map{assertKeyMatchRegex(it,"fasta",".*\\.(fasta|fa|fna)")}
        .map{it.fai!=null? it : it.plus(["fai":it.fasta+".fai"])}
        .map{it.dict!=null?it : it.plus(["dict":it.fasta.replaceAll("\\.(fasta|fa|fna)\$",".dict")])}
        .map{[file(it.fasta),file(it.fai),file(it.dict)]}
        .flatMap()
        .mix(Channel.of(fasta).map{it[1]})
        .mix(PREPARE_REFERENCE.out.fai.map{it[1]})
        .mix(PREPARE_REFERENCE.out.dict.map{it[1]})
        .map{fn->fn.toRealPath()}
        .unique()
        .collect()
        .map{[ref_hash,it]}



    bams_ch1 = bams_and_ref
            .map{if(!it.containsKey("bam")) throw new IllegalArgumentException("${it} : bam missing"); return it;}
            .map{if(!(it.bam.endsWith(".cram") || it.bam.endsWith(".bam"))) throw new IllegalArgumentException("${it}.bam should end with bam or cram"); return it;}
            .branch{
                /** sample undefined, need to find sample name */
                no_sample: !it.containsKey("sample") || it.sample.isEmpty() || it.sample.equals(".")
                /* sample defined */
                has_sample: true
            }
    
    /* get sample names for bam without id */
    SAMTOOLS_SAMPLES(
            all_references,
            bams_ch1.no_sample
                .map{file(it.bam)}
                .collect()
                .map{[[id:"bams"],it]}
            )
   versions = versions.mix(SAMTOOLS_SAMPLES.out.versions)

   fix_sample_name = SAMTOOLS_SAMPLES.out.samplesheet
            .map{meta,f->f}
            .splitCsv(sep:'\t',header:false)
            .map{[file(it[1]).toRealPath().toString(),it[0]]}
            .groupTuple()
            .map{bam,names->
                if(names.size()!=1) throw new IllegalArgumentException("Multiple samples for "+bam+":"+samples);
                return [bam,names[0]];
                }
            .join( bams_ch1.no_sample.map{[file(it.bam).toRealPath().toString(),it]})
            .map{bam,sample_name,meta->meta.plus("sample":sample_name)}

    bams_ch = fix_sample_name.mix(bams_ch1.has_sample)
			.map{it.bai?it: (it.bam.endsWith(".bam") ? it.plus(["bai":it.bam+".bai"]):  it.plus(["bai":it.bam+".crai"]))}
            .filter{assertKeyExistsAndNotEmpty(it,"sample")}
            .map{
		def meta = [id:it.sample];
		for(k in ["group","sex","father","mother","family","status","collection"])
			{
			if(!it.containsKey(k)) continue;
			def v = it.get(k);
			if(v.isEmpty()) continue;
			if(v.equals(".")) continue;
			meta=meta.plus([k:v]);
			}
		return [meta, file(it.bam),file(it.bai)];
		}
        
    /* check no duplicate samples */
    bams_ch.map{meta,bam,bai->meta.id}.unique().count()
        .combine(bams_ch.map{meta,bam,bai->meta.id}.count())
        .filter{c1,c2->c1!=c2}
        .view()
        .map{
            throw new IllegalArgumentException("Check the samplesheet. There is a duplicate sample name");
            }
    
    /** build pedigree from meta data */
    META_TO_PED(ref_hash, bams_ch.map{it[0]})
    versions = versions.mix(META_TO_PED.out.versions)


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
        .map{meta,beds->beds}
        .map{beds->beds instanceof List?beds:[beds]}
        .flatMap()
        .map{bed->[[id:bed.baseName],bed]}


    vcf_ch = Channel.empty()

    if(params.method.equalsIgnoreCase("gvcf")) {
        HAPLOTYPECALLER(
            [id:"hapcaller"],
            fasta,
            fai,
            dict,
            all_references,
            dbsnp,
            META_TO_PED.out.pedigree_gatk,
            beds_ch,
            bams_ch
            )
        versions = versions.mix(HAPLOTYPECALLER.out.versions)
        }
    else if(params.method.equalsIgnoreCase("bam2vcf")) {
        GATK_BAM2VCF(
            ref_hash,
            fasta,
            fai,
            dict,
            dbsnp,
            META_TO_PED.out.pedigree_gatk,
            all_references, //[meta, [ref files fa fai dict...]] all known reference
            beds_ch, // [meta,bed]
            bams_ch, // [meta,bam,bai]
            )
        versions = versions.mix(GATK_BAM2VCF.out.versions)
        vcf_ch = GATK_BAM2VCF.out.vcf
        }
    else if(params.method.equalsIgnoreCase("direct")) {
        HAPLOTYPECALLER_DIRECT(
            ref_hash,
            fasta,
            fai,
            dict,
            dbsnp,
            META_TO_PED.out.pedigree_gatk,
            beds_ch,
            bams_ch
            )
        versions = versions.mix(HAPLOTYPECALLER_DIRECT.out.versions)
        vcf_ch = HAPLOTYPECALLER_DIRECT.out.vcf
        }
    else
        {
        throw new IllegalArgumentException("undefined params.method=${params.method}")
        }
 
 /***************************************************
   *
   * BCFTOOLS STATS FROM VCFS
   *
   */
  BCFTOOLS_GUESS_PLOIDY(fasta, fai,vcf_ch)
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
    vcf_ch.map{[it[0],[it[1],it[2]]]}
    )
  versions = versions.mix(BCFTOOLS_GUESS_PLOIDY.out.versions)



	MULTIQC(
		ref_hash.plus("id":"hapcaller"),
		META_TO_PED.out.sample2collection,
		versions,
		multiqc
		)
    
}

runOnComplete(workflow)
