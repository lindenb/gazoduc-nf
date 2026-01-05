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

include { HAPLOTYPECALLER as HAPCALLER         }  from '../../../modules/gatk/hapcaller1'
include { BCFTOOLS_CONCAT                      }  from '../../../modules/bcftools/concat3'
include { COMBINE_GENOTYPE_GVCFS               }  from '../combinegenotypegvcfs'
include { makeKey                              }  from '../../../modules/utils/functions.nf'
include { flatMapByIndex                       }  from '../../../modules/utils/functions.nf'
include { FIND_GVCF_BLOCKS                     }  from '../../../modules/jvarkit/findgvcfblocks'
include { GLNEXUS_GENOTYPE                     }  from '../../../modules/glnexus/genotype'


workflow HAPLOTYPECALLER {
take:
    metadata
    fasta
    fai
    dict
    all_references// meta,fasta_files
    dbsnp // meta,vcf,tbi
    pedigree // meda,ped
    beds // meta,bed
    bams // meta,bam,bai
main:
    versions = Channel.empty()


	/* checl all beds have an unique ID */
	beds.map{meta,bed->[meta.id,bed]}
		.groupTuple()
		.filter{it[1].size()!=1}
		.map{throw new IllegalArgumentException("In HAPLOTYPECALLER bed should have a unique id"); return it;}

    HAPCALLER(
        fasta,
        fai,
        dict,
        all_references,
        bams.combine(beds)
            .map{meta1,bam,bai,meta2,bed->[meta1.plus(bed_id:meta2.id),bam,bai,bed]}
        )
    versions  = versions.mix(HAPCALLER.out.versions)

	/** group data by meta.bed_id */
    FIND_GVCF_BLOCKS(
        HAPCALLER.out.gvcf
            .map{meta,vcf,tbi,bed->[meta.bed_id,[vcf,tbi],bed]}
            .groupTuple()
            .map{bed_id,vcf_files,beds->[ [id:bed_id], vcf_files.flatten().sort(), beds.sort()[0]]}
        )
    versions  = versions.mix(FIND_GVCF_BLOCKS.out.versions)

    SPLIT_BED(FIND_GVCF_BLOCKS.out.bed)
    versions  = versions.mix(SPLIT_BED.out.versions)

    gvcfs_ch = SPLIT_BED.out.beds
		.map{meta,beds,_srcbed->[meta,(beds instanceof List?beds:[beds])]}
		.flatMap{flatMapByIndex(it,1)}
        .combine(HAPCALLER.out.gvcf) //join doesn't support duplicate keys, so use combine
		.filter{meta1,block_bed,meta2,vcf,tbi,srcbed->meta1.id==meta2.bed_id}
        .map{meta1,block_bed,meta2,vcf,tbi,srcbed->[meta1,vcf,tbi,block_bed]}
	
	
	vcf_out = Channel.empty()

	if(metadata.gvcf_merge_method==null) {
		log.warn("if(HAPLOTYPECALLER: gvcf_merge_method undefined")
		}

	if(metadata.gvcf_merge_method==null  || metadata.gvcf_merge_method.equalsIgnoreCase("combinegvcfs")) {
		COMBINE_GENOTYPE_GVCFS(
			metadata,
			fasta,
			fai,
			dict,
			dbsnp,
			gvcfs_ch
			)
		versions  = versions.mix(COMBINE_GENOTYPE_GVCFS.out.versions)
		vcf_out = COMBINE_GENOTYPE_GVCFS.out.vcf
		}
	else if(metadata.gvcf_merge_method.equalsIgnoreCase("glnexus")) {
		ch2 = gvcfs_ch
			.map{meta,gvcf,tbi,bed->
				[
				bed.toRealPath(),
				[gvcf,tbi]
				]}
			.groupTuple()
			.map{bed,vcf_files->[bed,vcf_files.flatten().sort()]}//need sort here to have the same key below
			.map{bed,vcf_files->[
				[id:makeKey([bed,vcf_files])],
				bed,
				vcf_files
				]}
			.multiMap{meta,bed_file,vcf_files->
				bed: [meta,bed_file]
				vcf: [meta,vcf_files]
				}
		GLNEXUS_GENOTYPE(
		        ch2.bed,
		       	[[id:"noconfig"],[]],
		        ch2.vcf
		        )
		versions  = versions.mix(GLNEXUS_GENOTYPE.out.versions)
		vcf_out = GLNEXUS_GENOTYPE.out.vcf
		}
	else
		{
		throw new IllegalArgumentException("unknown meta.gvcf_merge_method = ${metadata.gvcf_merge_method}.");
		}
    BCFTOOLS_CONCAT(
        vcf_out.map{meta,vcf,tbi,bed->[vcf,tbi]}//gvcf,tbi
			.flatMap()
            .collect()
            .map{files->[ [id:"gatkhapcaller"], files.sort() ]}
        )
    versions = versions.mix(BCFTOOLS_CONCAT.out.versions)
	vcf_out = BCFTOOLS_CONCAT.out.vcf
emit:
    versions
    vcf = vcf_out
}


process SPLIT_BED {
label "process_single"
label "array100"
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
stub:
"""
mkdir -p BEDS
touch BEDS/a.bed BEDS/b.bed BEDS/c.bed
touch versions.yml
"""
}
