
include {VCF_TO_BED         } from '../../modules/bcftools/vcf2bed'
include {BCFTOOLS_CONCAT    } from '../../modules/bcftools/concat'

workflow VQSR {
    take:
        meta
        fasta
        fai
        dict
        dbsnp
        recal_snps//val(meta),val(argument for recal snp)
        recal_indels//val(meta),val(argument for recal indels)
        vcfs
    main:
        versions = Channel.empty()

        BCF_TO_VCF(
            vcfs
                .map{[it[1],it[2]]}
                .collect()
                .map{[[id:"vqsr"],it]}
            )
        versions = versions.mix(BCF_TO_VCF.out.versions)

        COMPILE_MINIKIT(meta)
		versions = versions.mix(COMPILE_MINIKIT.out.versions)

        RECALIBRATE_SNP(
			fasta,
			fai,
			dict, 
			BCF_TO_VCF.out.vcf , 
			COMPILE_MINIKIT.out.jar,
			recal_snps
			)
		versions = versions.mix(RECALIBRATE_SNP.out.versions)

        RECALIBRATE_INDEL(
			fasta,
			fai,
			dict,
			BCF_TO_VCF.out.vcf , 
			COMPILE_MINIKIT.out.jar,
			recal_indels
			)
		versions = versions.mix(RECALIBRATE_INDEL.out.versions)



        
        vcfs = vcfs.map{[
            it[0].plus(id:it[1].toString().md5()),//give each vcf it's own id
            it[1],
            it[2]
            ]}
        
        
        VCF_TO_BED(vcfs)
        versions = versions.mix(VCF_TO_BED.out.versions)

        ch2 = VCF_TO_BED.out.output
            .splitCsv(sep:'\t',header:false,elem:1) //meta [ctg,start,end] vcf,vcfidx
            .map{[it[0], it[1][0]+":"+((it[1][1] as int)+1)+"-"+it[1][2], it[2], it[3]]} //meta, interval,vcf,vcfidx
        
		VCF_TO_INTERVALS(ch2)
		versions = versions.mix(VCF_TO_INTERVALS.out.versions)
        
       
        ch3 = VCF_TO_INTERVALS.out.bed
            .combine(RECALIBRATE_SNP.out.vcf.map{[it[1]/*snp*/,it[2]/*snpidx*/,it[3]/* tranche */]})
            .combine(RECALIBRATE_INDEL.out.vcf.map{[it[1]/*index*/,it[2]/*indelidx*/,it[3]/* tranche */]})
            

        APPLY_RECALIBRATION_SNP(
            fasta,
			fai,
			dict,
            ch3
            )
		versions = versions.mix(APPLY_RECALIBRATION_SNP.out.versions)

		APPLY_RECALIBRATION_INDEL(
            fasta,
			fai,
			dict,
            APPLY_RECALIBRATION_SNP.out.vcf
            )
		versions = versions.mix(APPLY_RECALIBRATION_INDEL.out.versions)
		
        BCFTOOLS_CONCAT(
            APPLY_RECALIBRATION_INDEL.out.vcf
                .map{[[id:"vqsr"],[it[1],it[2]]]}
                .groupTuple()
                .map{[it[0],it[1].flatten()]},
            [[id:"nobed"],[]]
            )
        versions = versions.mix(BCFTOOLS_CONCAT.out.versions)
    emit:
        versions
        vcf = BCFTOOLS_CONCAT.out.vcf
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
	def prefix = task.ext.prefix?:"reformat"+(meta.id?:"")
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
}


process VCF_TO_INTERVALS {
tag "${interval} ${vcf.name}"
label "process_single"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta),val(interval),path(vcf),path(idx)
output:
	tuple val(meta),path("*.bed"),path(vcf),path(idx),emit:bed
	path("versions.yml"),emit:versions
script:
	def distance="10Mb"
	def min_distance=100
    def prefix = interval.md5().substring(0,7)
"""
hostname 1>&2
mkdir -p TMP
bcftools view -G "${vcf}" "${interval}" |\\
	jvarkit -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP vcf2intervals \
		--bed \
		--distance "${distance}" \
		--min-distance "${min_distance}" > ${prefix}.bed

cat << END_VERSIONS > versions.yml
"${task.process}":
	java: todo
	jvarkit: todo
END_VERSIONS
"""
}


process COMPILE_MINIKIT {
label "process_short"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	val(meta)
output:
	tuple val(meta),path("minikit.jar"),emit:jar
	path("versions.yml"),emit:versions
script:
"""

mkdir -p TMP/TMP
cp -v "${moduleDir}/Minikit.java" TMP/Minikit.java

cat << EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: Minikit
EOF

javac -d TMP -sourcepath TMP TMP/Minikit.java
jar cfm minikit.jar TMP/tmp.mf -C TMP .


cat << END_VERSIONS > versions.yml
"${task.process}":
	java: todo
END_VERSIONS
"""
}


process RECALIBRATE_SNP {
tag "${vcf.name}"
label "process_medium"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta),path(vcf),path(vcfidx)
	tuple val(meta4),path(minikit)
	tuple val(meta5),val(recal_snps)
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),path("*.tranches.txt"),emit:vcf
	path("*.plot.R"),emit:R
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:vcf.name.md5()+".snp.recal"
	def atts = task.ext.atts?:"ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum"
	if(recal_snps.trim().isEmpty()) throw new IllegalArgumentException("missing  args in ${task.process}");
"""
hostname 1>&2
set -x
mkdir -p TMP


bcftools view --type snps -G "${vcf}" |\\
	java  -Djava.io.tmpdir=TMP  -XX:-UsePerfData -jar ${minikit} --an '${atts}' > TMP/args.list
test -s TMP/args.list

cat << EOF >> TMP/args.list
${recal_snps}
EOF

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData" VariantRecalibrator  \\
	-R "${fasta}" \\
	-V "${vcf}" \\
	--arguments_file TMP/args.list \\
	-mode SNP \\
	-O "${prefix}.vcf.gz" \\
	--tranches-file "${prefix}.tranches.txt" \\
    --dont-run-rscript \\
	--rscript-file  "${prefix}.plot.R"

cat << END_VERSIONS > versions.yml
"${task.process}":
	java: todo
	jvarkit: todo
END_VERSIONS
"""
}


process RECALIBRATE_INDEL {
tag "${vcf.name}"
label "process_medium"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta),path(vcf),path(vcfidx)
	tuple val(meta4),path(minikit)
	tuple val(meta5),val(recal_indels)
output:
	tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi"),path("*.tranches.txt"),emit:vcf
	path("*.plot.R"),emit:R
	path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:vcf.name.md5()+".indel.recal"
    def atts = task.ext.atts?:"ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum"

"""
hostname 1>&2
mkdir -p TMP
set -x

bcftools view --type indels -G "${vcf}" |\\
	java  -Djava.io.tmpdir=TMP   -XX:-UsePerfData -jar ${minikit} --an ${atts} > TMP/args.list

test -s TMP/args.list

cat << EOF  >> TMP/args.list
${recal_indels}
EOF

# 20200630 add AS https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
# finalement AS ne semble pas marcher avec indel

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantRecalibrator  \\
	-R "${fasta}" \\
	-V "${vcf}" \\
	--arguments_file TMP/args.list \\
	-mode INDEL \\
	-O "${prefix}.recal.vcf.gz" \\
	--tranches-file "inde${prefix}.tranches.txt" \\
    --dont-run-rscript \\
	--rscript-file ${prefix}.plot.R

cat << END_VERSIONS > versions.yml
"${task.process}":
	java: todo
	jvarkit: todo
END_VERSIONS
"""
}

process APPLY_RECALIBRATION_SNP {
tag "${vcf.name} ${bed.name}"
label "process_subgle"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta),
        path(bed),
        path(vcf),path(vcfidx),
        path(recal_snp_vcf),path(recal_snp_vcf_idx),path(recal_snp_tranches),
        path(recal_indel_vcf),path(recal_indel_vcf_idx),path(recal_indel_tranches)
output:
	tuple val(meta),
        path(bed),
		path("*.vcf.gz"),path("*.vcf.gz.tbi"),
		path(recal_indel_vcf),path(recal_indel_vcf_idx),path(recal_indel_tranches),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:(bed.name+vcf.name).md5().substring(0,7)+".snp"
    def level = task.ext.level?:99.0
"""
hostname 1>&2
mkdir -p TMP


if ${vcf.name.endsWith(".bcf")} 
then
    bcftools view --threads ${task.cpus} --regions-file "${bed}"  -O z -o TMP/jeter.vcf.gz ${vcf}
    bcftools index --threads ${task.cpus} -f -t  TMP/jeter.vcf.gz
fi

gatk  --java-options "-Xmx${task.memory.giga}g  -XX:-UsePerfData -Djava.io.tmpdir=TMP" ApplyVQSR  \\
        -R "${fasta}" \\
        -V ${vcf.name.endsWith(".bcf")?"TMP/jeter.vcf.gz":"${vcf}"} \\
        -mode SNP \\
        -L "${bed}" \\
        --truth-sensitivity-filter-level ${level} \\
        --recal-file "${recal_snp_vcf}" \\
        --tranches-file "${recal_snp_tranches}" \\
        -O "TMP/jeter2.vcf.gz"

mv TMP/jeter2.vcf.gz ${prefix}.vcf.gz
mv TMP/jeter2.vcf.gz.tbi ${prefix}.vcf.gz.tbi

cat << END_VERSIONS > versions.yml
"${task.process}":
	java: todo
	bcftools: todo
END_VERSIONS
"""
}


process APPLY_RECALIBRATION_INDEL {
tag "${vcf.name} ${bed.name}"
label "process_subgle"
afterScript 'rm -rf  TMP'
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
	tuple val(meta1),path(fasta)
	tuple val(meta2),path(fai)
	tuple val(meta3),path(dict)
	tuple val(meta ),path(bed),path(vcf),path(vcf_idx),
		path(recal_indel_vcf),path(recal_indel_vcf_idx),path(recal_indel_tranches)
output:
	tuple val(meta ),path("*.vcf.gz"),path("*.vcf.gz.tbi"),emit:vcf
    path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:(bed.name+vcf.name).md5().substring(0,7)+".indel"
    def level = task.ext.level?:99.0
"""
hostname 1>&2
mkdir -p TMP

gatk  --java-options "-Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP" ApplyVQSR  \\
	-R "${fasta}" \\
	-V "${vcf}"  \\
	-mode INDEL \\
	-L "${bed}" \\
	--truth-sensitivity-filter-level ${level} \\
	--recal-file "${recal_indel_vcf}" \\
	--tranches-file "${recal_indel_tranches}" \\
	-O "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz "${prefix}.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${prefix}.vcf.gz.tbi"


cat << END_VERSIONS > versions.yml
"${task.process}":
	java: todo
	jvarkit: todo
END_VERSIONS
"""
}
