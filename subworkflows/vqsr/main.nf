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
include { BCFTOOLS_CONCAT                  } from '../../modules/bcftools/concat3'
include { verify                           } from '../../modules/utils/functions.nf'
include { isBlank                          } from '../../modules/utils/functions.nf'
include { flatMapByIndex                   } from '../../modules/utils/functions.nf'
include { VARIANT_RECALIBRATOR             } from '../../modules/gatk/variantrecalibrator'
include { APPLY_VQSR  as APPLY_VQSR_SNPS   } from '../../modules/gatk/applyvqsr'
include { APPLY_VQSR  as APPLY_VQSR_INDELS } from '../../modules/gatk/applyvqsr'
include { SPLIT_N_VARIANTS                 } from '../../modules/jvarkit/splitnvariants'
include { BCFTOOLS_INDEX                   } from '../../modules/bcftools/index'

workflow VQSR {
    take:
        metadata
        fasta
        fai
        dict
        dbsnp
        vcfs
    main:
        versions = Channel.empty()
		multiqc = Channel.empty()

		/** extract concatenated variant of VCF input without genotypes */
        BCF_TO_VCF(
            vcfs.map{meta,vcf,tbi->[vcf,tbi]}
				.flatMap()
                .collect()
                .map{files->[metadata,files.sort()]}
            )
        versions = versions.mix(BCF_TO_VCF.out.versions)

        GUESS_PARAMETERS(
			BCF_TO_VCF.out.vcf.combine(Channel.of("snps","indels"))
				.map{meta,vcf,idx,t->[meta.plus(type:t),vcf,idx]}
			)
		versions = versions.mix(GUESS_PARAMETERS.out.versions)

		/** transfert meta2.type to meta */
		ch1 = BCF_TO_VCF.out.vcf
			.combine(GUESS_PARAMETERS.out.arguments)
			.filter{meta,vcf,tbi,meta2,args->meta.id==meta2.id} //paranoid
			.multiMap{meta,vcf,tbi,meta2,args->
				arguments : [meta2,args]
				vcf :  [meta.plus(type:meta2.type),vcf,tbi]
			}

        VARIANT_RECALIBRATOR(
			fasta,
			fai,
			dict,
			ch1.arguments,
			ch1.vcf 
			)
		versions = versions.mix(VARIANT_RECALIBRATOR.out.versions)

		SPLIT_N_VARIANTS(
			[[id:"nobed"],[]],
			vcfs
			)
		versions = versions.mix(SPLIT_N_VARIANTS.out.versions)
		
		vcfs = SPLIT_N_VARIANTS.out.vcf.flatMap(row->flatMapByIndex(row,1))
			.combine(SPLIT_N_VARIANTS.out.tbi.flatMap(row->flatMapByIndex(row,1)))
			.filter{meta1,vcf,meta2,tbi->meta1.id==meta2.id && "${vcf.name}.tbi" == tbi.name}
			.map{meta1,vcf,meta2,tbi->[meta1.plus(id:vcf.name.md5()),vcf,tbi]}
		
		ch1 = vcfs.combine(VARIANT_RECALIBRATOR.out.vcf.filter{meta,vcf,tbi,tranches->meta.type=="snps"})
			.multiMap{meta1,vcf,tbi,meta2,recal,recalidx,tranches->
				recal : [meta2,recal,recalidx,tranches]
				vcf :  [meta1,vcf,tbi]
			}


		APPLY_VQSR_SNPS(
			fasta,
			fai,
			dict,
			ch1.recal,
			ch1.vcf
			)
		versions = versions.mix(APPLY_VQSR_SNPS.out.versions)

		
		ch1 = APPLY_VQSR_SNPS.out.vcf
			.combine(VARIANT_RECALIBRATOR.out.vcf.filter{meta,vcf,tbi,tranches->meta.type=="indels"})
			.multiMap{meta1,vcf,tbi,meta2,recal,recalidx,tranches->
				recal : [meta2,recal,recalidx,tranches]
				vcf :  [meta1,vcf,tbi]
			}

		APPLY_VQSR_INDELS(
			fasta,
			fai,
			dict,
			ch1.recal,
			ch1.vcf
			)
		versions = versions.mix(APPLY_VQSR_INDELS.out.versions)

		/* for index by bcftools to get index metadata */
		BCFTOOLS_INDEX(APPLY_VQSR_INDELS.out.vcf.map{meta,vcf,tbi->[meta,vcf]})
		versions = versions.mix(BCFTOOLS_INDEX.out.versions)

    emit:
        versions
		multiqc
        vcf = BCFTOOLS_INDEX.out.vcf
}


process BCF_TO_VCF {
tag "${meta.id?:""}"
label "process_short"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript 'rm -rf  TMP'
input:
	tuple val(meta),path("VCFS/*")
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
	tuple val(meta),path("chroms.txt"),emit:contigs
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"reformat${meta.id?:""}"
"""
hostname 1>&2

mkdir -p TMP
find VCFS \\( -name "*.vcf.gz" -o -name "*.bcf" \\) > TMP/jeter.list
test -s TMP/jeter.list

bcftools concat --threads ${task.cpus} --drop-genotypes -a -O z -o TMP/${prefix}.vcf.gz --file-list TMP/jeter.list
bcftools index  --threads ${task.cpus} --force --tbi TMP/${prefix}.vcf.gz 

mv TMP/${prefix}.* ./
	
bcftools index --stats ${prefix}.vcf.gz | cut -f 1 > chroms.txt

cat << END_VERSIONS > versions.yml
"${task.process}":
	bcftools: todo
END_VERSIONS
"""
stub:
	def prefix = task.ext.prefix?:"reformat${meta.id?:""}"
 """
touch chroms.txt ${prefix}.vcf.gz ${prefix}.vcf.gz.tbi versions.yml
"""
}





process GUESS_PARAMETERS {
tag "${meta.type}"
label "process_medium"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta),path(vcf),path(vcfidx)
output:
	tuple val(meta),path("*.list"),emit:arguments
	path("versions.yml"),emit:versions
script:
	def jvm = task.ext.jvm?:"-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData"
	def type = meta.type
	verify(!isBlank(type),"type is blank")
	def attributes = task.ext.attributes?:""
	verify(!isBlank(attributes),"atts is blank. Should be like 'ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum' ")
	def prefix = task.ext.prefix?:"${meta.id}.${type}"

	def other_recal_args = task.ext.recal_args?:""
	verify(!isBlank(other_recal_args),"For ${type} missing ext.recal_args in ${task.process}")

"""
mkdir -p TMP
cp -v "${moduleDir}/Minikit.java" TMP/Minikit.java
javac -d TMP -sourcepath TMP TMP/Minikit.java


bcftools view --type ${type} -G "${vcf}" |\\
	java  ${jvm} -cp TMP Minikit --an '${attributes}' > TMP/args.list

test -s TMP/args.list

cat << EOF >> TMP/args.list
${other_recal_args}
EOF

mv TMP/args.list "${prefix}.list"

touch versions.yml
"""
stub:
	def type = meta.type
	verify(!isBlank(type),"type is blank")
	def prefix = task.ext.prefix?:"${meta.id}.${type}"
"""
touch versions.yml "${prefix}.list"
"""
}


