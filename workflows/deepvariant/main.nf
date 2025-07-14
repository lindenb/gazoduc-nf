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

include {DEEPVARIANT_TRIOS} from '../../subworkflows/deepvariant/trios/main.nf'

include {runOnComplete; dumpParams} from '../../modules/utils/functions.nf'


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
	def bed =      [ refhash, (params.bed==null?[]:file(params.bed))]
	def pedigree = [ refhash, (params.pedigree==null?[]:file(params.pedigree))]
	
	versions = Channel.empty()

	if(params.bed==null) {
		SCATTER_TO_BED(refhash,fasta,fai,dict)
		versions = versions.mix(SCATTER_TO_BED.out.versions)
		bed = SCATTER_TO_BED.out.bed
		}


	BEDCLUSTER(fasta,fai,dict,bed)
	versions = versions.mix(BEDCLUSTER.out.versions)

	bed = BEDCLUSTER.out.bed.flatMap{
		def L=[]
		for(f in it[1]) {
			L.add([[id:f.name],f]);
			}
		return L;
		}
	
/*
	
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
		.map{[[id:it.sample],file(it.bam),file(it.bai)]}

	ch1 = bams_ch.combine(bed).multiFlatMap{
		bam: [it[0],it[1],it[3]],
		bed:[it[3],it[4]
		}


	DEEPVARIANT_TRIOS(
		[id:"trios"],
		fasta,
		fai,
		dict,
		pedigree,
		ch1.bed,
		ch1.bam
		)
	versions = versions.mix(TRIOS.out.versions)

	BCFTOOLS_CONCAT(
		TRIOS.out.bed
			.map{[it[1],it[2]]}
			.collect()
			.flatten()
			.map{[[id:"concat"],it]}}
		)
*/	
		
	}

runOnComplete(workflow)


