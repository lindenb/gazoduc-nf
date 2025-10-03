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

include { HAPLOTYPECALLER as HAPCALLER         }  from '../../../modules/gatk/hapcaller1'
include { BCFTOOLS_CONCAT                      }  from '../../../modules/bcftools/concat'
include { COMBINE_GENOTYPE_GVCFS               }  from '../combinegenotypegvcfs'
include { makeKey                              }  from '../../../modules/utils/functions.nf'
include {FIND_GVCF_BLOCKS                      }  from '../../../modules/jvarkit/findgvcfblocks'


workflow HAPLOTYPECALLER {
take:
    meta
    fasta
    fai
    dict
    all_references// meta,fasta_files
    dbsnp // meta,vcf,tbi
    beds // meta,bed
    bams // meta,bam,bai
main:
    versions = Channel.empty()

    HAPCALLER(
        fasta,
        fai,
        dict,
        all_references,
        bams.combine(beds)
            .map{meta1,bam,bai,meta2,bed->[meta1,bam,bai,bed]}
        )
    versions  = versions.mix(HAPCALLER.out.versions)


    FIND_GVCF_BLOCKS(
        HAPCALLER.out.gvcf
            .map{meta,vcf,tbi,bed->[bed.toRealPath(),[vcf,tbi]]}
            .groupTuple()
            .map{bed,vcf_files->[ [id:"${bed.name}"], vcf_files.flatten().sort(), bed]}
        )
    versions  = versions.mix(FIND_GVCF_BLOCKS.out.versions)

    SPLIT_BED(FIND_GVCF_BLOCKS.out.bed)
    versions  = versions.mix(SPLIT_BED.out.versions)

    gvcfs_ch = SPLIT_BED.out.beds
		.map{meta,beds,srcbed->[meta,(beds instanceof List?beds:[beds]),srcbed]}
		.flatMap{meta,beds,srcbed->{
			def L=[];
			for(int i=0;i< beds.size();i++) {
				L.add([srcbed,beds[i]]);
				}
			return L;
			}}
        .combine(HAPCALLER.out.gvcf) // I tried to use 'join' but it doesn't work / Nasty bug ? So combine+filter
        .filter{srcbed1,bed,meta,vcf,tbi,srcbed2->srcbed1.toRealPath().equals(srcbed2.toRealPath())}
        .map{srcbed1,bed,meta,vcf,tbi,srcbed2->[meta,vcf,tbi,bed]}


    COMBINE_GENOTYPE_GVCFS(
            meta,
            fasta,
            fai,
            dict,
	        dbsnp,
            gvcfs_ch
            )

    versions  = versions.mix(COMBINE_GENOTYPE_GVCFS.out.versions)


    BCFTOOLS_CONCAT(
        COMBINE_GENOTYPE_GVCFS.out.vcf
            .map{meta,vcf,tbi,bed->[vcf,tbi]}//gvcf,tbi
             .collect()
             .map{[[id:"gatkhapcaller"],it]},
        [[:],[]] //bed
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

emit:
    versions
    vcf = BCFTOOLS_CONCAT.out.vcf
}


process SPLIT_BED {
tag "${meta.id} ${bed.name} ${src_bed}"
afterScript "rm -rf TMP"
conda  "${moduleDir}/../../../conda/bioinfo.01.yml"
input:
        tuple val(meta),path(bed),path(src_bed)
output:
        tuple val(meta),path("BEDS/*.bed",arity:'1..*'),path(src_bed),emit:beds
	path("versions.yml"),emit:versions
script:
        def f = meta.id
        def size=task.ext.size?:500000;
"""
mkdir -p BEDS TMP

sort -T TMP -t '\t' -k1,1 -k2,2n "${bed}" |\\
awk -F '\t' 'BEGIN{
	N=1;
	T=0.0;
	f=sprintf("BEDS/${f}.%d.N${size}.bed",N);
	}
	{
	print \$0 >> f;
	T+=int(\$3)-int(\$2);
	if(T>=${size}) {
		close(f);
		T=0.0;
		N++;
		f=sprintf("BEDS/${f}.%d.N${size}.bed",N);
		}
	}'

# remove oberlaps
find BEDS -type f -name "*.bed" | while read B
do
	sort -T TMP -t '\t' -k1,1 -k2,2n "\${B}" |\\
		bedtools merge > TMP/jeter.bed && mv TMP/jeter.bed "\${B}"
done

touch versions.yml
"""
}
