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

include {SMOOVE_SV                     } from '../../subworkflows/smoove/population' 
include {SMOOVE_ANNOTATE               } from '../../modules/smoove/annotate' 
include {JVARKIT_VCF_SET_DICTIONARY    } from '../../modules/jvarkit/vcfsetdict' 
include {runOnComplete                 } from '../../modules/utils/functions.nf'
include {k1_signature                  } from '../../modules/utils/k1.nf'



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
            ucsc_name: (params.ucsc_name?:"undefined"),
			ensembl_name: (params.ensembl_name?:"undefined")
            ]
        def fasta  = [ref_hash, file(params.fasta) ]
        def fai    = [ref_hash, file(params.fai) ]
        def dict   = [ref_hash, file(params.dict) ]
		def gff3    = [ref_hash, file(params.gff3) , file(params.gff3+".tbi")]
		def exclude    = [[id:"no_exclude"], [] ]


	if(params.exclude_bed==null) {
		GET_EXCLUDE(fasta,fai,dict)
		versions  = versions.mix(GET_EXCLUDE.out.versions)
		exclude = GET_EXCLUDE.out.bed
		}
	else {
		exclude    = [ref_hash, file(params.exclude_bed) ]
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
		.map{[[id:it.sample],file(it.bam),file(it.bai)]}
           
	if(params.cases!=null) {
		//all are controls by default
		bams = bams.map{[it[0].id,it[0].plus("status":"control"),it[1],it[2]]}

		// set case for children in trio or if status is affected
		samples_status = Channel.fromPath(params.cases)
			.splitText(S->S.trim())
			.map{[it,"case"]}
			
		
		// join bams, 
		bams = bams.join(samples_status, failOnMismatch:false,remainder:true )
			.filter{it[1]!=null}
			.map{
				def hash=it[1];
				if(it[4]!=null && it[4].equals("case")) {
					hash = hash.plus("status":"case");
					}
				return [hash,it[2],it[3]]
				}
			.unique()
		}

	SMOOVE_SV(
		[id:"smoove"],
		fasta,
		fai,
		dict,
		exclude,
		bams
		)
	versions  = versions.mix(SMOOVE_SV.out.versions)

	vcf = SMOOVE_SV.out.vcf
	if((params.with_annotation as boolean)==true) {
		SMOOVE_ANNOTATE(fasta,fai,dict,gff3,vcf)
		versions  = versions.mix(SMOOVE_ANNOTATE.out.versions)
		vcf = SMOOVE_ANNOTATE.out.vcf
		}
	
	
	JVARKIT_VCF_SET_DICTIONARY(dict,vcf)
	versions  = versions.mix(JVARKIT_VCF_SET_DICTIONARY.out.versions)
	vcf = JVARKIT_VCF_SET_DICTIONARY.out.vcf

	}



process GET_EXCLUDE {
tag "${fasta.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
    tuple val(meta3),path(dict)
output:
	tuple val(meta1),path("*.bed"),emit:bed /* exclude2.bed otherwise collision with bed from scatter_bed */
    path("versions.yml"),emit:versions
script:
 def k1 = k1_signature()
 def prefix = task.ext.prefix?:fasta.name+".exclude"
"""
hostname 1>&2

mkdir -p TMP

#cf https://github.com/brentp/smoove

cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter1.tsv
1:${k1.hg19}\thttps://github.com/hall-lab/speedseq/blob/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed
1:${k1.hg38}\thttps://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed
EOF


awk -F '\t' '{printf("%s:%s\\n",\$1,\$2);}' '${fai}' |\\
	sed 's/^chr//' |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter2.tsv

join -t '\t' -1 1 -2 1 -o '1.2' TMP/jeter1.tsv TMP/jeter2.tsv |\\
	sort | uniq > TMP/jeter.url

test -s TMP/jeter.url

url=`cat  TMP/jeter.url`

curl -L "\${url}" |\\
	cut -f 1,2,3|\\
	jvarkit bedrenamechr -R "${fasta}" --column 1 --convert SKIP > TMP/exclude2.bed 

awk -F '\t' '(!(\$1 ~ /^(chr)?[0-9XY]+\$/)) {printf("%s\t0\t%s\\n",\$1,\$2);}' '${fai}' >> TMP/exclude2.bed

sort -T TMP -t '\t' -k1,1 -k2,2n TMP/exclude2.bed > ${prefix}.bed

test -s ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
    url: \${url}
END_VERSIONS
"""
}


runOnComplete(workflow)
