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

*/
include {MULTIQC                                  } from '../../modules/multiqc'
include {COMPILE_VERSIONS                         } from '../../modules/versions/main.nf'
include {VCF_STATS                                } from '../../subworkflows/vcfstats'
include {PIHAT                                    } from '../../subworkflows/pihat'
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
		versions = Channel.empty()
		multiqc = Channel.empty()

		def ref_hash = [
            id: file(params.fasta).simpleName,
            ucsc_name: (params.ucsc_name?:"undefined")
            ]
    def fasta  = [ref_hash, file(params.fasta) ]
    def fai    = [ref_hash, file(params.fai) ]
    def dict   = [ref_hash, file(params.dict) ]
		def gtf    = [ref_hash, file(params.gtf), file(params.gtf+".tbi") ]
		def gff3    = [ref_hash, file(params.gff3), file(params.gff3+".tbi") ]
    def bed =  [ref_hash, []]
    def pedigree = [ref_hash, []]
    if(params.bed!=null) {
        bed =  [ref_hash, file(params.bed)];
    }
    if(params.pedigree!=null) {
        pedigree =  [ref_hash, file(params.pedigree)];
    }

    ch1 = Channel.empty()
    if(params.samplesheet.endsWith(".list")) {
      ch1 = Channel.fromPath(params.samplesheet)
          .splitText()
          .map{it.trim()}
          .map{[vcf:it]}
          }
    else if(params.samplesheet.endsWith(".vcf.gz") || params.samplesheet.endsWith(".bcf")) {
        ch1 = Channel.of(params.samplesheet).map{[vcf:it]}
      }
    else {
        ch1 = Channel.fromPath(params.samplesheet).splitCsv(header:true,sep:',')
        }

    vcfs = ch1
			.map{assertKeyMatchRegex(it,"vcf","^\\S+\\.(bcf|vcf\\.gz)\$")}
		.map{
				if(it.containsKey("index")) return it;
				if(it.vcf.endsWith(".vcf.gz")) return it.plus(index : it.vcf+".tbi");
				return it.plus(index:it.vcf+".csi");
			}
      .map{
				if(it.containsKey("id")) return it;
				return it.plus(id:it.vcf.toString().md5());
			}
			.map{[
        [id:it.id],
        file(it.vcf),
        file(it.index)
      ]}

    if(params.onekgenome==null) {
      ch1 = Channel.empty()
      }
   else if(params.onekgenome.endsWith(".list")) {
      ch1 = Channel.fromPath(params.onekgenome)
          .splitText()
          .map{it.trim()}
          .filter{it.length()>0}
          .map{[vcf:it]}
      }
    else {
        ch1 = Channel.fromPath(params.onekgenome).splitCsv(header:true,sep:',')
        }

    onekgenome = ch1
        .map{assertKeyMatchRegex(it,"vcf","^\\S+\\.(bcf|vcf\\.gz)\$")}
        .map{
          if(it.containsKey("index")) return it;
          if(it.vcf.endsWith(".vcf.gz")) return it.plus(index : it.vcf+".tbi");
          return it.plus(index:it.vcf+".csi");
        }
        .map{
          if(it.containsKey("id")) return it;
          return it.plus(id:it.vcf.toString().md5());
        }
        .map{[
          [id:it.id],
          file(it.vcf),
          file(it.index)
        ]}


    PIHAT(
       	ref_hash,
        fasta,
        fai,
        dict,
        onekgenome,
        vcfs
    )
	versions = versions.mix(PIHAT.out.versions)


	 VCF_STATS(
			ref_hash,
			fasta,
			fai,
			dict,
			Channel.of(gtf),
			Channel.of(gff3),
			pedigree,
			[[id:"nogroup2sample"],[]],
			bed,
			vcfs
			)
	versions = versions.mix(VCF_STATS.out.versions)
	multiqc = versions.mix(VCF_STATS.out.multiqc)

	COMPILE_VERSIONS(versions.collect())
	multiqc = multiqc.mix(COMPILE_VERSIONS.out.multiqc)

  MULTIQC(multiqc.collect().map{[[id:"vcfstats"],it]})
	}
runOnComplete(workflow)

