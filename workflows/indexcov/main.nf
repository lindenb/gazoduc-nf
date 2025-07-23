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


include {INDEXCOV                   } from '../../subworkflows/indexcov/simple'
include {JVARKIT_BAM_RENAME_CONTIGS } from '../../modules/jvarkit/bamrenamechr'
include {runOnComplete;dumpParams   } from '../../modules/utils/functions.nf'
include {COMPILE_VERSIONS           } from '../../modules/versions'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


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
	versions = Channel.empty()
 	def ref_hash = [
            id: file(params.fasta).simpleName,
            ucsc_name: (params.ucsc_name?:"undefined")
            ]
	def fasta  = [ref_hash, file(params.fasta) ]
	def fai    = [ref_hash, file(params.fai) ]
	def dict   = [ref_hash, file(params.dict) ]
	def pedigree =  [ref_hash,[] ]
	if(params.pedigree!=null) {
		pedigree = [ref_hash, file(params.pedigree)];
		}

  	bams = Channel.fromPath(params.samplesheet)
			.splitCsv(header:true,sep:',')
			.map{assertKeyMatchRegex(it,"sample","^[A-Za-z_0-9\\.\\-]+\$")}
			.map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
			.map{
				if(it.containsKey("bai")) return it;
				if(it.bam.endsWith(".cram")) return it.plus(bai : it.bam+".crai");
				return it.plus(bai:it.bam+".bai");
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
			.branch {
				ok_ref: it.fasta.equals(params.fasta)
				bad_ref: true
				}
			//[[],]}
            //.groupTuple() // group by fasta
	
	
	JVARKIT_BAM_RENAME_CONTIGS(
		dict,
		[[id:"nobed"],[]],
		bams.bad_ref.map{[[id:it.sample],file(it.bam),file(it.bai),file(it.fasta),file(it.fai),file(it.dict)]}
		)
	versions =versions.mix(JVARKIT_BAM_RENAME_CONTIGS.out.versions)

	def collate_size = (params.batch_size as int)

	INDEXCOV(
		ref_hash.plus([batch_size: collate_size, id:"indexcov"]),
		fasta,
		fai,
		dict,
		pedigree,
		bams.ok_ref.map{[[id:it.sample],file(it.bam),file(it.bai)]}
			.mix(JVARKIT_BAM_RENAME_CONTIGS.out.bam.map{[it[0].plus([force_rebuild:false]),it[1],it[2]]})
		)
	versions = versions.mix(INDEXCOV.out.versions)

	COMPILE_VERSIONS(versions.collect())
	}


runOnComplete(workflow)
