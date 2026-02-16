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
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include { COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include { VCF_STATS                                } from '../../subworkflows/vcfstats'
include { PIHAT                                    } from '../../subworkflows/pihat'
include { runOnComplete; dumpParams                } from '../../modules/utils/functions.nf'
include { VCF_INPUT                                } from '../../subworkflows/nf/vcf_input'
include { VCF_INPUT  as VCF_INPUT_1KG              } from '../../subworkflows/nf/vcf_input'
include { PREPARE_ONE_REFERENCE                    } from '../../subworkflows/samtools/prepare.one.ref'
include { MULTIQC as MQC1                          } from '../../modules/multiqc'
include { MULTIQC as MQC2                          } from '../../modules/multiqc'
include { JVARKIT_MULTIQCPOSTPROC                  } from '../../modules/jvarkit/multiqcpostproc'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}

workflow {
		versions = Channel.empty()
		multiqc = Channel.empty()

		def metadata = [id:"vcfstats"]
          
    if(params.vcf==null) {
		log.warn("undefined --vcf")
		exit -1
		}
	if(params.fasta==null) {
		log.warn("undefined --fasta")
		exit -1
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




		def gtf    = [metadata, file(params.gtf), file(params.gtf+".tbi") ]
		def gff3    = [metadata, file(params.gff3), file(params.gff3+".tbi") ]
    def pedigree = [metadata, []]
    def sample2pop = [metadata, []]

    if(params.bed!=null) {
        bed =  [metadata, file(params.bed)];
      }
    else
      {
      bed = PREPARE_ONE_REFERENCE.out.scatter_bed
      }

    if(params.pedigree!=null) {
        pedigree =  [metadata, file(params.pedigree)];
      }

    if(params.sample2pop!=null) {
        sample2pop  =  [metadata, file(params.sample2pop)];
    }

	/***************************************************
	*
	* VCF input
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
	
	vcfs = VCF_INPUT.out.vcf
	    .map{meta,vcf,tbi->[metadata,vcf,tbi]}


  if(params.with_pihat==true) {
      
      VCF_INPUT_1KG(metadata.plus([
        path: params.onekgenome,
        arg_name: "onekgenome",
        require_index : true,
        required: true,
        unique : false
        ]))
      versions = versions.mix(VCF_INPUT_1KG.out.versions)
      
      PIHAT(
          metadata,
          fasta,
          fai,
          dict,
          [[id:"nosample2pop"],[]],
          [[id:"noexcludesamples"],[]],
          [[id:"noexcludebed"],[]],
          VCF_INPUT_1KG.out.vcf,
          vcfs
      )
    versions = versions.mix(PIHAT.out.versions)
    }

	 VCF_STATS(
			metadata.plus(
        with_bcftools_stats:true,
        ad_ratio : true,
        gatk_denovo : true
        ),
			PREPARE_ONE_REFERENCE.out.fasta,
			PREPARE_ONE_REFERENCE.out.fai,
			PREPARE_ONE_REFERENCE.out.dict,
			Channel.of(gtf),
			Channel.of(gff3),
			pedigree,
			sample2pop,
			bed,
			vcfs
			)
	versions = versions.mix(VCF_STATS.out.versions)
	multiqc = multiqc.mix(VCF_STATS.out.multiqc)

	COMPILE_VERSIONS(versions.collect())
	multiqc = multiqc.map{meta,f->f}.mix(COMPILE_VERSIONS.out.multiqc)


  MQC1(
    [[id:"no_mqc_config"],[]],
    multiqc.collect().map{[[id:"vcfstats"],it.sort()]}
    )

	JVARKIT_MULTIQCPOSTPROC(
		sample2pop,
		[[id:"nocustom"],[]],
		MQC1.out.datadir
		)
	

	MQC2(
		[[id:"no_mqc_config"],[]],
		JVARKIT_MULTIQCPOSTPROC.out.json.map{_meta,f->f}.collect().map{files->[metadata.plus(id:"perpop"),files.sort()]}
		)

	}
runOnComplete(workflow)

