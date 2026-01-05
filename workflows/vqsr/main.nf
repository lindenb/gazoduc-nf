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

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
  // log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
//validateParameters()

// Print summary of supplied parameters
//log.info paramsSummaryLog(workflow)

include {VQSR                                     } from '../../subworkflows/vqsr'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include {VCF_STATS                                } from '../../subworkflows/vcfstats'
include {runOnComplete; dumpParams                } from '../../modules/utils/functions.nf'



if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
		versions = Channel.empty()
		multiqc = Channel.empty()

		def ref_hash = [
            id: file(params.fasta).simpleName,
            ucsc_name: (params.ucsc_name?:"undefined")
            ]
        def fasta  = [ref_hash, file(params.fasta) ]
        def fai    = [ref_hash, file(params.fai) ]
        def dict   = [ref_hash, file(params.dict) ]
		def dbsnp  = [ref_hash, file(params.dbsnp) , file(params.dbsnp+".tbi") ] 
		def gtf    = [ref_hash, file(params.gtf), file(params.gtf+".tbi") ]
		def gff3    = [ref_hash, file(params.gff3), file(params.gff3+".tbi") ]

		def recal_snps = [ref_hash, params.recal_snps_args]
		def recal_indels = [ref_hash,params.recal_indels_args]

		if(params.vcf.endsWith(".list")) {
			in_vcf_ch = Channel.fromPath(params.vcf).
				splitText().
				map{it.trim()}.
				map{file(it)}
			}
		else
			{
			in_vcf_ch = Channel.fromPath(params.vcf)
			}
		in_vcf_ch = in_vcf_ch
			.map{[
				[id:it.toString().md5()],
				file(it), 
				file(it.name.endsWith(".bcf")?""+it+".csi":""+it+".tbi")
				]}
	
		
	
	VQSR(
		ref_hash,
		fasta,
		fai,
		dict,
		dbsnp,
		recal_snps,
		recal_indels,
		in_vcf_ch
		)
	versions = versions.mix(VQSR.out.versions)

	
	 VCF_STATS(
			ref_hash,
			fasta,
			fai,
			dict,
			Channel.of(gtf),
			Channel.of(gff3),
			[[id:"noped"],[]],
			[[id:"nogroup2sample"],[]],
			[[id:"nobed"],[]],
			VQSR.out.vcf
			)
	versions = versions.mix(VCF_STATS.out.versions)
	multiqc = versions.mix(VCF_STATS.out.multiqc)

	COMPILE_VERSIONS(versions.collect())
	multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

    MULTIQC(multiqc.collect().map{[[id:"vqsr"],it]})

	}
runOnComplete(workflow)
