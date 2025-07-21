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

include {dumpParams;runOnComplete  } from '../../modules/utils/functions.nf'
include { SCATTER_TO_BED           } from '../../subworkflows/gatk/scatterintervals2bed'
include { MANTA_MULTI              } from '../../modules/manta/multi'
include { MANTA_CONVERT_INVERSION  } from '../../modules/manta/convert.inversion'
include { TRUVARI_COLLAPSE         } from '../../modules/truvari/collapse'
include { ANNOTATE_SV              } from '../../subworkflows/annotation/sv'
include { JVARKIT_VCF_SET_DICTIONARY } from '../../modules/jvarkit/vcfsetdict'


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

    def ref_hash = [
        id: file(params.fasta).simpleName,
        ucsc_name: (params.ucsc_name?:"undefined"),
        ensembl_name: (params.ensembl_name?:"undefined"),
        ]
    def fasta  = [ref_hash, file(params.fasta) ]
    def fai    = [ref_hash, file(params.fai) ]
    def dict   = [ref_hash, file(params.dict) ]
    def gtf   = [ref_hash, file(params.gtf), file(params.gtf+".tbi") ]

    versions = Channel.empty()


    ch1 = Channel.fromPath(params.samplesheet)
        .splitCsv(header:true,sep:',')
        .map{assertKeyMatchRegex(it,"sample","^[A-Za-z_0-9\\.\\-]+\$")}
        .map{assertKeyMatchRegex(it,"bam","^\\S+\\.(bam|cram)\$")}
        .map{
            if(it.containsKey("batch")) return it;
            return it.plus(batch:"batch_"+it.sample);
        }
        .map{assertKeyMatchRegex(it,"batch","^[A-Za-z_0-9\\.\\-]+\$")}
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
        .map{assertKeyMatchRegex(it,"bai","^\\S+\\.(bai|crai)\$")}
        .map{[ it.batch, it.sample,file(it.bam),file(it.bai),file(it.fasta),file(it.fai),file(it.dict)]}
        .groupTuple()
        .map{
            java.util.HashSet<String> set = new  java.util.HashSet<String>();
            set.addAll(it[4]);
            if(set.size()!=1) throw new IllegalArgumentException("multiple fasta file for ${it} for the same batch.")
            return it;
        }
        .multiMap {
            fasta: [ref_hash,it[4][0]]
            fai: [ref_hash,it[5][0]]
            dict: [ref_hash,it[6][0]]
            bams: [ref_hash.plus(id:it[0]), it[2].plus(it[3]) ]
        }
   
                
    if(params.bed==null) {
        SCATTER_TO_BED(ref_hash,fasta,fai,dict)
        versions = versions.mix(SCATTER_TO_BED.out.versions)
        bed = SCATTER_TO_BED.out.bed
        }
	else {
        bed =  [ref_hash, file(params.bed) ]
        }

    MANTA_MULTI(
        ch1.fasta,
        ch1.fai,
        ch1.dict,
        bed,
        ch1.bams
        )
    versions = versions.mix(MANTA_MULTI.out.versions.first())
    
    JVARKIT_VCF_SET_DICTIONARY(
        dict,
        MANTA_MULTI.out.diploidSV
        )
    versions = versions.mix(JVARKIT_VCF_SET_DICTIONARY.out.versions.first())


    MANTA_CONVERT_INVERSION(
        fasta,
        fai,
        dict,
        JVARKIT_VCF_SET_DICTIONARY.out.vcf
        )
    versions = versions.mix(MANTA_CONVERT_INVERSION.out.versions.first())

    TRUVARI_COLLAPSE(
        fasta,
        fai,
        dict,
        MANTA_CONVERT_INVERSION.out.vcf
            .map{[it[1],it[2]]}
            .collect()
            .map{[[id:"truvari"],it.flattenn()]}
        )
    versions = versions.mix(TRUVARI_COLLAPSE.out.versions.first())
	
    ANNOTATE_SV(
        ref_hash,
        fasta,
        fai,
        dict,
        gtf,
        TRUVARI_COLLAPSE.out.vcf
        )
    versions = versions.mix(ANNOTATE_SV.out.versions.first())
    }


runOnComplete(workflow)
