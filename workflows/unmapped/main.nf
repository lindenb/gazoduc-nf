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
include {UNMAPPED                                 } from '../../subworkflows/unmapped'
include {runOnComplete;dumpParams                 } from '../../modules/utils/functions.nf'
include { SAMTOOLS_STATS                          } from '../../modules/samtools/stats'
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'



Map assertKeyExists(final Map hash,final String key) {
    if(!hash.containsKey(key)) throw new IllegalArgumentException("no key ${key}'in ${hash}");
    return hash;
}

Map assertKeyExistsAndNotEmpty(final Map hash,final String key) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(value.isEmpty()) throw new IllegalArgumentException("empty ${key}'in ${hash}");
    return hash;
}

Map assertKeyMatchRegex(final Map hash,final String key,final String regex) {
    assertKeyExists(hash,key);
    def value = hash.get(key);
    if(!value.matches(regex)) throw new IllegalArgumentException(" ${key}'in ${hash} doesn't match regex '${regex}'.");
    return hash;
}




workflow {
	def versions = Channel.empty()
	def to_multiqc  = Channel.empty()
	def meta = [id:"contaminations"]
	def acn_file = [meta,file("${moduleDir}/contaminants.txt")];
	if(params.acn_file!=null) {
		acn_file = [[id:file(params.acn_file).baseName],file(params.acn_file)]
	}


	bams_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true,sep:',')
        .map{assertKeyMatchRegex(it,"sample","^[A-Za-z_0-9\\.\\-]+\$")}
        .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
		.map{
            if(it.containsKey("bai")) return it;
            if(it.bam.endsWith(".cram")) return it.plus(bai : it.bam+".crai");
            return it.plus(bai:it.bam+".bai");
        	}
		.map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
		.map{
				if(it.containsKey("fasta")) return it;
				if(params.fasta==null || params.fasta.trim().isEmpty()) throw new IllegalArgumentException("undefined --fasta");
				return it.plus(fasta:params.fasta);
			}
			.map{
				if(it.containsKey("fai")) return it;
				return it.plus(fai:it.fasta+".fai");
			}
			.map{[[id:it.sample], file(it.bam), file(it.bai), file(it.fasta), file(it.fai)]}


	ch2 = bams_ch.multiMap{
		fasta: [[id:"ref"],it[3]]
		fai  : [[id:"ref"],it[4]]
		bam :  [it[0],it[1],it[2]]
		}
	

	SAMTOOLS_STATS(
		ch2.fasta,
		ch2.fai,
	    [[id:"NOBED"],[]],
		ch2.bam
	    )
	versions = versions.mix(SAMTOOLS_STATS.out.versions)
	to_multiqc = to_multiqc.mix(SAMTOOLS_STATS.out.stats.map{it[1]})

	UNMAPPED(
		meta,
		acn_file,
		bams_ch
		)
	versions = versions.mix(UNMAPPED.out.versions)
	to_multiqc = to_multiqc.mix(UNMAPPED.out.multiqc)

 	COMPILE_VERSIONS(versions.collect())
	to_multiqc = to_multiqc.mix(COMPILE_VERSIONS.out.multiqc)

    MULTIQC(to_multiqc.collect().map{[[id:"contaminations"],it]})
	}


runOnComplete(workflow);
