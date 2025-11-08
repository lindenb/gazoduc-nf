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


include {dumpParams                 } from '../../modules/utils/functions.nf'
include {runOnComplete              } from '../../modules/utils/functions.nf'
include {WGSIM                      } from '../../subworkflows/samtools/wgsim'
include {MULTIQC                    } from '../../modules/multiqc'
include {COMPILE_VERSIONS           } from '../../modules/versions'





workflow {

	if( params.help ) {
    	dumpParams(params);
    	exit 0
		} 
	else {
		dumpParams(params);
		}

	if(params.fasta==null) {
		log.err("undefined --fasta");
		exit -1
	}


   def metadata = [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
      ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ metadata, file(params.fasta)]
	versions = Channel.empty()
	multiqc_ch = Channel.empty()

	
	
	WGSIM(
		metadata.plus([
			n_samples: (params.n_samples as int),
			random_seed : (params.random_seed as long)
			]),
		Channel.of(fasta)
		)
	versions = versions.mix(WGSIM.out.versions)

	MAKE_SAMPLESHEET(
		WGSIM.out.fastq
		.map{meta,R1,R2->[
			meta.id,
			"${params.outdir}/FASTQ/${params.prefix}${R1.name}",
			"${params.outdir}/FASTQ/${params.prefix}${R2.name}",
			meta.sex?:"",
			meta.father?:"",
			meta.mother?:"",
			meta.status?:"",
			meta.collection?:""
			]}
		.map{it.join(",")}
		.collect()
		)
	versions = versions.mix(MAKE_SAMPLESHEET.out.versions)
	
	
    COMPILE_VERSIONS(versions.collect().map{it.sort()})
    multiqc_ch = multiqc_ch.mix(COMPILE_VERSIONS.out.multiqc.map{[[id:"versions"],it]})
    // in case of problem multiqc_ch.filter{!(it instanceof List) || it.size()!=2}.view{"### FIX ME ${it} MULTIQC"}
    MULTIQC(
	[[id:"no_mqc_config"],[]],
	multiqc_ch.map{it[1]}.collect().map{[[id:"wgsim"],it]}
	)
    }


runOnComplete(workflow);


process MAKE_SAMPLESHEET {
tag "N=${L.size()}"
input:
	val(L)
output:
	path("samplesheet.csv"),emit:samplesheet
	path("versions.yml"),emit:versions
script:
"""

echo 'sample,fastq_1,fastq_2,sex,father,mother,status,collection' > jeter.csv
cat << EOF >> jeter.csv
${L.join("\n")}
EOF

mv jeter.csv samplesheet.csv

touch versions.yml
"""
}
