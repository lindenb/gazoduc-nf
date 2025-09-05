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
LIABILITY, WH
ETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


ORIGINAL snakemake workflow by Raphael Blanchet PhD.

*/
include {runOnComplete; dumpParams                } from '../../modules/utils/functions.nf'

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

def hash_ref= [
		id: file(params.fasta).baseName,
		name: file(params.fasta).baseName,
    ucsc_name: (params.ucsc_name?:"undefined")
		]
	def fasta = [ hash_ref, file(params.fasta)]
	def fai   = [ hash_ref, file(params.fai)]
	def dict  = [ hash_ref, file(params.dict)]



		splitCsv(sep:'\t',header:true).
		map{[ [:], it.sample, it.R1, it.R2]}

  	samplesheet_ch =Channel.fromPath(params.samplesheet).
			.map{assertKeyMatchRegex(it,"sample","^[A-Za-z][A-Za-z0-9_\\-\\.]*[A-Za-z0-9]\$")}
			.map{assertKeyMatchRegex(it,"fastq_1","^\\S+\\.(fastq|fq|ora|fastq\\.gz|fq\\.gz)\$")}
      .map{
          if(it.containsKey("fastq_2")) {
            if(!(it.fastq_2.equals(".") || it.fastq_2.isEmpty())) {
              // fastq cannot end with '.ora'
              assertKeyMatchRegex(it,"fastq_2","^\\S+\\.(fastq|fq|fastq\\.gz|fq\\.gz)\$")
              }
            return it;
            }
          }
			

	
	}
