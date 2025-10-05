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
include {PREPARE_REFERENCE          } from '../../subworkflows/samtools/prepare.ref'
include {WGSIM                      } from '../../modules/samtools/wgsim'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}





workflow {
  def hash_ref= [
      id: file(params.fasta).baseName,
      name: file(params.fasta).baseName,
      ucsc_name: (params.ucsc_name?:"undefined")
      ]
	def fasta = [ hash_ref, file(params.fasta)]
	
	PREPARE_REFERENCE(hash_ref,fasta)
	versions = versions.mix(PREPARE_REFERENCE.out.versions)
	fai = PREPARE_REFERENCE.out.fai

	
	def n_samples= (params.n_samples as int)
	ch1 = Channel.of(0..<n_samples).map{n->[id:"S"+(n+1)]}
	
	WGSIM(fasta,fai,ch1)
	versions = versions.mix(WGSIM.out.versions)
	
	MAKE_SAMPLESHEET(WGSIM.out.fastq.map{meta,R1,R2->[
			meta.id,
			"${params.outdir}/FASTQ/${params.prefix}${R1.name}",
			"${params.outdir}/FASTQ/${params.prefix}${R2.name}",
			]}
		.map{it.join(",")}
		)
	versions = versions.mix(MAKE_SAMPLESHEET.out.versions)
	
    }


runOnComplete(workflow);


process MAKE_SAMPLESHEET {
input:
	val(L)
output:
	path("samplesheet.csv"),emit:samplesheet
	path("versions.yml"),emit:versions
script:
"""

echo 'sample,fastq_1,fastq_2' > jeter.csv
cat << EOF >> jeter.csv
${L.join("\n")}
EOF

mv jeter.csv samplesheet.csv

touch versions.yml
"""
}
