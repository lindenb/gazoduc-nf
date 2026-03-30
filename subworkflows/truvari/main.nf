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
include { TRUVARI_COLLAPSE       } from '../../modules/truvari/collapse'
include { TRUVARI_ANNO_CHUNCKS   } from '../../modules/truvari/anno.chuncks'
include { verify                 } from '../../modules/utils/functions'
include { isBlank                } from '../../modules/utils/functions'
include { parseBoolean           } from '../../modules/utils/functions'
include { makeKey                } from '../../modules/utils/functions'
include { flatMapByIndex         } from '../../modules/utils/functions'
include { BCFTOOLS_CONCAT        } from '../../modules/bcftools/concat3'


workflow TRUVARI {
     take:
		metadata
		fasta
		fai
		dict
        vcfs /* [meta,vcf] */
     main:
	 	versions = Channel.empty()
		multiqc = Channel.empty()




		BCFTOOLS_SPLIT_BY_TYPE(
			fai,
			vcfs
			)
		versions = versions.mix(BCFTOOLS_SPLIT_BY_TYPE.out.versions)
	
		typed_vcf_ch = BCFTOOLS_SPLIT_BY_TYPE.out.symbolic.map{meta,vcf,tbi->[meta.plus(variant_type:"symbolic"),vcf,tbi]}
			.mix(BCFTOOLS_SPLIT_BY_TYPE.out.bnd.map{meta,vcf,tbi->[meta.plus(variant_type:"bnd"),vcf,tbi]})
			.mix(BCFTOOLS_SPLIT_BY_TYPE.out.atgc.map{meta,vcf,tbi->[meta.plus(variant_type:"atgc"),vcf,tbi]})
			.map{meta,vcf,tbi->[meta.plus(id:"${meta.id}.${meta.variant_type}"),vcf,tbi]}
	
		BCFTOOLS_MERGE_FOR_TRUVARI(
			fasta,
			fai,
			typed_vcf_ch
				.map{meta,vcf,tbi->[
					[
					id : meta.variant_type,
					variant_type:meta.variant_type
					] ,[vcf,tbi]]
					}
				.groupTuple()
				.map{meta,array->[meta,array.flatten().sort()]}
			)
		versions = versions.mix(BCFTOOLS_MERGE_FOR_TRUVARI.out.versions)


		/*
		If you find spans in the counts.bed with a huge number of SVs, these are prime candidates for exclusion to speed up truvari collapse. To exclude them, you first subset to the regions of interest and then use bedtools to create a bed file that will skip them.
		*/
		if(metadata.with_anno_chunks==null || parseBoolean(metadata.with_anno_chunks)) {
			TRUVARI_ANNO_CHUNCKS(
				[[id:"nobed"],[]],
				BCFTOOLS_MERGE_FOR_TRUVARI.out.vcf.map{meta,vcf,tbi->[meta,vcf]}
				)
			versions = versions.mix(TRUVARI_ANNO_CHUNCKS.out.versions)
			}
		
		ch1  = BCFTOOLS_MERGE_FOR_TRUVARI.out.bed
			.flatMap{row->flatMapByIndex(row,1)}
			.combine(BCFTOOLS_MERGE_FOR_TRUVARI.out.vcf)
			.filter{meta1,_bed,meta2,_vcf,_tbi->meta1.id==meta2.id}
			.multiMap{meta1,bed,meta2,vcf,tbi->
				bed: [meta1,bed]
				vcf : [meta2.plus(id:makeKey([meta2,vcf,bed])),vcf,tbi]
			}

		TRUVARI_COLLAPSE(
			fasta,
			fai,
			ch1.bed,
			ch1.vcf
			)
		versions = versions.mix(TRUVARI_COLLAPSE.out.versions)
		
		BCFTOOLS_CONCAT(
			TRUVARI_COLLAPSE.out.vcf
				.flatMap{_meta,vcf,tbi->[vcf,tbi]}
				.collect()
				.map{files->[metadata,files.sort()]}
			)
		versions = versions.mix(BCFTOOLS_CONCAT.out.versions)

	emit:
		vcf = BCFTOOLS_CONCAT.out.vcf
		versions
		multiqc
	}


process BCFTOOLS_SPLIT_BY_TYPE {
tag "${meta.id}"
label "process_single"
label "array100"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript 'rm -rf  TMP'
input:
		tuple val(meta ),path(fai)
        tuple val(meta ),path(vcf),path(tbi)
output:
        tuple val(meta),path("*.symbolic.bcf"),path("*.symbolic.bcf.csi"),emit:symbolic
		tuple val(meta),path("*.atgc.bcf"),path("*.atgc.bcf.csi"),emit:atgc
		tuple val(meta),path("*.bnd.bcf"),path("*.bnd.bcf.csi"),emit:bnd
        path("versions.yml"),emit:versions
script:
	def args1 = task.ext.args1?:""
	def prefix = task.ext.prefix?:"${meta.id}"
	def awk_expr = task.ext.awk_expr?:"(\$1 ~ /^(chr)?[0-9XY]+\$/)"
"""
mkdir -p TMP

awk -F '\t' '${awk_expr} {printf("%s\t0\t%s\\n",\$1,\$2);}' '${fai}' |\\
	sort -T TMP -k1,1V -k2,2n > TMP/jeter.bed

test -s TMP/jeter.bed

bcftools view \\
	${args1} \\
	--targets-file TMP/jeter.bed \\
	--threads ${task.cpus} \\
	-O v -o TMP/jeter.vcf \\
	"${vcf}" 

# see https://github.com/ACEnglish/truvari/wiki/collapse#symbolic-variants
# wee need to help bcftools to manage combine with ID
 
 awk -F '\t' '/^#/ {print;next;} (\$5 ~ /^</) {print;}' TMP/jeter.vcf |\\
	bcftools annotate --set-id '%CHROM:%POS:%END:%SVTYPE' -O b -o TMP/symbolic.bcf

awk -F '\t' '/^#/ {print;next;} (!(\$5 ~ /^</) && !(\$5 ~ /^[ATGCN]+\$/ ) ) {print;}' TMP/jeter.vcf |\\
	bcftools annotate --set-id '%CHROM:%POS:%FIRST_ALT:%SVTYPE' -O b -o TMP/bnd.bcf

awk -F '\t' '/^#/ {print;next;} (!(\$5 ~ /^</) && (\$5 ~ /^[ATGCN]+\$/ )) {print;}' TMP/jeter.vcf |\\
	bcftools view -O b -o TMP/atgc.bcf

for F in symbolic atgc bnd
do
	bcftools index -f TMP/\${F}.bcf
	mv TMP/\${F}.bcf ${prefix}.\${F}.bcf
	mv TMP/\${F}.bcf.csi ${prefix}.\${F}.bcf.csi
done


cat << __EOF__ > versions.yml
${task.process}:
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
__EOF__
"""

stub:
def prefix = task.ext.prefix?:"${meta.id}"
"""
touch ${prefix}.symbolic.bcf ${prefix}.symbolic.bcf.csi  \\
	${prefix}.atgc.bcf ${prefix}.atgc.bcf.csi  \\
	${prefix}.bnd.bcf ${prefix}.atgc.bnd.csi  \\
	versions.yml
"""
}

process BCFTOOLS_MERGE_FOR_TRUVARI {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript 'rm -rf  TMP'
input:
        tuple val(meta1),path(fasta)
		tuple val(meta2),path(fai)
        tuple val(meta ),path("VCFS/*")
output:
        tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
		tuple val(meta),path("contigs*.bed", arity: '0..*'),emit:bed
        path("versions.yml"),emit:versions
script:
		def variant_type = meta.variant_type
		verify(!isBlank(variant_type),"variant_type is blank")
		def args1 = task.ext.args1?:""
		def args2 = task.ext.args2?:"--filter-logic '+'"
		def args3 = task.ext.args3?:""
		def merge_mode = task.ext.merge_mode?:(variant_type=="atgc"?"none":"id")
		def prefix = task.ext.prefix?:"${meta.id}.merge"
"""
hostname 1>&2
mkdir -p TMP
set -x
find VCFS/ -type l \\( -name "*.vcf.gz" -o -name "*.bcf" \\) |\\
		while read F
		do
				bcftools query -l "\${F}" | head -n 1 | tr "\\n" "\t" >> TMP/jeter.txt
				echo "\${F}" >> TMP/jeter.txt
		done

sort -T TMP -t '\t' -k1,1 TMP/jeter.txt | cut -f 2 > TMP/jeter.list


bcftools merge \\
		--threads ${task.cpus} \\
		${args1} \\
		${args2} \\
		--force-samples \\
		--file-list TMP/jeter.list \\
		-m ${merge_mode} \\
		-O b \\
		-o TMP/jeter2.bcf
	
mv TMP/jeter2.bcf TMP/jeter.bcf

#  Truvari: Cannot compare multi-allelic records. Please split
if ${merge_mode!="id"}
then
bcftools norm \\
	--threads ${task.cpus} \\
	--fasta-ref ${fasta} \\
	--multiallelics -any \\
	-O u \\
	-o TMP/jeter2.bcf \\
	TMP/jeter.bcf

	mv TMP/jeter2.bcf TMP/jeter.bcf
fi


if ${!isBlank(args3)}
then
	
	# optional filter ?
	bcftools view --threads ${task.cpus} ${args3}  -O u -o TMP/jeter2.bcf TMP/jeter.bcf
	
	mv TMP/jeter2.bcf TMP/jeter.bcf
fi



# bug DRAGEN
bcftools annotate --threads ${task.cpus} --force -x 'FORMAT/SR' -O z -o TMP/jeter2.vcf.gz TMP/jeter.bcf
bcftools index --threads ${task.cpus} -f -t TMP/jeter2.vcf.gz

bcftools index -s TMP/jeter2.vcf.gz |\\
	awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' |\\
	split -a 9 --lines=1 --additional-suffix=.bed - TMP/contigs.${meta.id}

mv -v TMP/contigs.*bed ./ || true

mv TMP/jeter2.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter2.vcf.gz.tbi ${prefix}.vcf.gz.tbi




cat << __EOF__ > versions.yml
${task.process}:
    bcftools: "\$(bcftools version | awk '(NR==1) {print \$NF;}')"
__EOF__
"""

stub:
	def prefix = task.ext.prefix?:"${meta.id}.merge"
"""
touch ${prefix}.bcf ${prefix}.bcf.csi  versions.yml
"""
}

/**

When considering what threshold you would like to use, just looking at the 4th column may not be sufficient as the 'sub-chunks' may be smaller and will therefore run faster. There also is no guarantee that a region with high SV density will be slow. For example, if all SVs in a chunk could collapse, it would only take O(N - 1) comparisons to complete the collapse. The problems arise when the SVs have few redundant SVs and therefore requires O(N**2) comparisons.

*/
process ANNO_CHUNCKS_COMPLEMENT {
tag "${meta.id?:""}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript 'rm -rf  TMP'
input:
		tuple val(meta1),path(fai)
        tuple val(meta ),path(bed)
output:
        tuple val(meta),path("*.bed"),emit:bed
        path("versions.yml"),emit:versions
script:
	def max_sv = task.ext.max_sv?:30000
	def prefix = task.ext.prefix?:"${meta.id}.select${max_sv}"
"""
mkdir -p TMP

awk '\$4 >= ${max_sv}' "${bed}" |\\
	sort -t '\t' -k1,1 -k2,2n > TMP/to_exclude.bed


if test ! -s TMP/to_exclude.bed
then
	tail -1 '${fai}' | awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' > TMP/to_exclude.bed
fi

awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' '${fai}' |\\
	sort -t '\t' -k1,1 -k2,2n > TMP/genome.bed

bedtools complement -g TMP/genome.bed -i TMP/to_exclude.bed > TMP/jeter.bed
mv TMP/jeter.bed ${prefix}.bed

touch versions.yml
"""

stub:
	def max_sv = task.ext.max_sv?:30000
	def prefix = task.ext.prefix?:"${meta.id}.select${max_sv}"
"""
touch versions.yml ${prefix}.bed
"""
}