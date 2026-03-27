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
include { verify                 } from '../../modules/utils/functions'
include { isBlank                } from '../../modules/utils/functions'
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
		def regex_contig = metadata.regex_contig?:"(chr)?[0-9XY]+"




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

		TRUVARI_COLLAPSE(
			fasta,
			fai,
			[[id:"nobed"],[]],
			BCFTOOLS_MERGE_FOR_TRUVARI.out.vcf
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