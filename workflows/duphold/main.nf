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
include {MULTIQC                     } from '../../modules/multiqc'
include {COMPILE_VERSIONS            } from '../../modules/versions/main.nf'
include {runOnComplete; dumpParams   } from '../../modules/utils/functions.nf'
include {DUPHOLD as RUN_DUPHOLD      } from '../../subworkflows/duphold'


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


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}



workflow {
	def refhash=[
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
		ucsc_name :( params.ucsc_name?:"undefined"),
		ensembl_name : (params.ensembl_name?:"undefined")
		]

	def fasta =    [ refhash, file(params.fasta) ]
	def fai   =    [ refhash, file(params.fai)  ]
	def dict  =    [ refhash, file(params.dict) ]
	def gtf  =     [ refhash, file(params.gtf) ]
	def vcf   =    [ refhash, file(params.vcf), file(params.vcf + (params.vcf.endsWith(".vcf.gz")?".tbi":".csi"))]
	def pedigree = [ refhash, []]
	def snv_vcf  = [ refhash, [], []]


	if(params.snv_vcf!=null) {
		snv_vcf  =  [ refhash, file(params.snv_vcf), file(params.snv_vcf+(params.snv_vcf.endsWith(".vcf.gz")?".tbi":".csi"))]
		}
	if(params.pedigree!=null) {
		pedigree  =  [ refhash, file(params.pedigree) ]
		}
	versions = Channel.empty()
	multiqc = Channel.empty()


  

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
				if(it.containsKey("batch")) return it;
				return it.plus(batch:it.sample);
			}
		.map{
				if(it.containsKey("fasta")) return it;
				return it.plus(fasta:params.fasta);
			}
			.map{
				if(it.containsKey("fai")) return it;
				return it.plus(fai:it.fasta+".fai");
			}
			.map{
				if(it.containsKey("dict")) return it;
				return it.plus(dict: it.fasta.replaceAll("\\.(fasta|fa|fna)\$",".dict"));
			}
		.map{[
			[id: it.sample, fasta:it.fasta.md5()],
			file(it.bam),
			file(it.bai),
			file(it.fasta),
			file(it.fai),
			file(it.dict)
			]}
	
	
	RUN_DUPHOLD(
		[id:"duphold"],
		fasta,
		fai,
		dict,
		bams_ch,
		Channel.of(vcf),
		Channel.of(snv_vcf)
		)
	versions = versions.mix(RUN_DUPHOLD.out.versions)
	/*

	COMPILE_VERSIONS(versions.collect())
	multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)
	
	MULTIQC(
		[[id:"no_mqc_config"],[]],
		multiqc.collect().map{[[id:"bcftools"],it]})
	*/
	}

runOnComplete(workflow)



