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
include {BCFTOOLS_CONCAT} from '../../subworkflows/bcftools/bcftools.concat.02.nf'
include {moduleLoad; escapeXml; dumpParams; runOnComplete} from '../../modules/utils/functions.nf'
include {JVARKIT_VCF_TO_INTERVALS_01} from '../../subworkflows/jvarkit/jvarkit.vcf2intervals.nf'

if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}


workflow {
	VARIANT_VQSR_38(
		params.genomeId,
		params.vcf,
		file(params.recal_snps),
		file(params.recal_indels),
		file(params.bed)
		)
	}


runOnComplete(workflow)

workflow VARIANT_VQSR_38 {
	take:
		genomeId
		vcf
		recal_snps
		recal_indels
		output_bed
	main:


		vcf2bed_ch = BCF_TO_VCF(vcf)

		compile_ch= COMPILE_MINIKIT()


		vcf2rgn_ch = JVARKIT_VCF_TO_INTERVALS_01([distance:"10Mb",min_distance:"100"],vcf,output_bed)

		rgn_vcf_ch = vcf2rgn_ch.bed.splitCsv(sep:'\t',header:false).
			map{ [it[0]+":"+((it[1] as int)+1)+"-"+it[2], it[3] ]}


		recal_snp_ch = RECALIBRATE_SNP( genomeId, vcf2bed_ch.output,compile_ch.output, vcf2rgn_ch.vcf2bed, recal_snps)

		recal_indel_ch = RECALIBRATE_INDEL(genomeId, vcf2bed_ch.output ,compile_ch.output, vcf2rgn_ch.vcf2bed, recal_indels)


		

		apply_snp_ch = APPLY_RECALIBRATION_SNP(genomeId , recal_snp_ch.vcf, recal_snp_ch.tranches, rgn_vcf_ch)
		apply_indel_ch = APPLY_RECALIBRATION_INDEL(genomeId , recal_indel_ch.vcf, recal_indel_ch.tranches, apply_snp_ch.output)


		concat_ch = BCFTOOLS_CONCAT([method:"all"],apply_indel_ch.output.map{[vcf:it.toString()]},file("NO_FILE"))
		
/*
	emit:
		version = version_ch
		vcf = gather_ch.vcf
		index = gather_ch.index
*/
	}


process BCF_TO_VCF {
tag "${file(vcf).name}"
afterScript 'rm -rf  TMP'
input:
	val(vcf)
output:
	path("reformat.vcf.gz"),emit:output
	path("reformat.vcf.gz.tbi"),emit:index
	path("chroms.txt"),emit:contig
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
mkdir -p TMP

if ${vcf.endsWith(".list")} ; then

	bcftools concat --threads ${task.cpus} --drop-genotypes -a -O z -o TMP/reformat.vcf.gz --file-list "${vcf}"
	bcftools index  --threads ${task.cpus}  --tbi TMP/reformat.vcf.gz 
	mv TMP/reformat.* ./
	
elif ${vcf.endsWith(".bcf")} ; then

	bcftools view  --threads ${task.cpus} --drop-genotypes -O z -o TMP/reformat.vcf.gz "${vcf}"
	bcftools index  --threads ${task.cpus} --tbi TMP/reformat.vcf.gz 
	mv TMP/reformat.* ./

else

	ln -s "${vcf}" reformat.vcf.gz
	ln -s "${vcf}.tbi" reformat.vcf.gz.tbi
fi

bcftools index --stats reformat.vcf.gz | cut -f 1 > chroms.txt

"""
}


process COMPILE_MINIKIT {
executor "local"
afterScript "rm -rf TMP"
output:
	path("minikit.jar"),emit:output
script:
"""
${moduleLoad("java/1.8")}

mkdir -p TMP/TMP
cp -v "${moduleDir}/Minikit.java" TMP/Minikit.java

cat << EOF > TMP/tmp.mf
Manifest-Version: 1.0
Main-Class: Minikit
EOF

javac -d TMP -sourcepath TMP TMP/Minikit.java
jar cfm minikit.jar TMP/tmp.mf -C TMP .
"""
}


process RECALIBRATE_SNP {
tag "file(vcf).name"
afterScript 'rm -rf  TMP'
input:
	val(genomeId)
	val(vcf)
	path(minikit)
	path(bed)
	path(recal_snps)
output:
	path("${params.prefix}snp.recal.vcf.gz"),emit:vcf
	path("${params.prefix}snp.tranches.txt"),emit:tranches
	path("${params.prefix}snp.recal.vcf.gz.tbi")
	path("${params.prefix}snp.tranches.pdf")
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("java/1.8 r bcftools")}

mkdir -p TMP

test -s "${recal_snps}"

bcftools view --regions-file "${bed}" --type snps -G "${vcf}" | java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum > TMP/args.list
test -s TMP/args.list


cat ${recal_snps} >> TMP/args.list

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
	-T VariantRecalibrator \\
	-R "${reference}" \\
	-L "${bed}" \\
	--input "${vcf}" \\
	--arg_file TMP/args.list \\
	-mode SNP \\
	--recal_file "TMP/snp.recal.vcf.gz" \\
	--tranches_file "TMP/snp.tranches.txt" \\
	-rscriptFile  "TMP/snp.plot.R"

mv -v "TMP/snp.recal.vcf.gz" "${params.prefix}snp.recal.vcf.gz"
mv -v "TMP/snp.recal.vcf.gz.tbi" "${params.prefix}snp.recal.vcf.gz.tbi"
mv -v "TMP/snp.tranches.txt" "${params.prefix}snp.tranches.txt"
mv -v "TMP/snp.tranches.txt.pdf" "${params.prefix}snp.tranches.pdf"

"""
}



process RECALIBRATE_INDEL {
tag "${file(vcf).name}"
afterScript 'rm -rf  TMP'
input:
	val(genomeId)
	val(vcf)
	path(minikit)
	path(bed)
	path(recal_indels)
output:
	path("${params.prefix}indel.recal.vcf.gz"),emit:vcf
	path("${params.prefix}indel.tranches.txt"),emit:tranches
	path("${params.prefix}indel.recal.vcf.gz.tbi")
script:
	def reference = params.genomes[genomeId].fasta

"""
hostname 1>&2
${moduleLoad("java/1.8 r bcftools")}
mkdir -p TMP

test -s "${recal_indels}"

bcftools view --regions-file "${bed}" --type indels -G "${vcf}" | java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum > TMP/args.list
test -s TMP/args.list


cat ${recal_indels} >> TMP/args.list

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
	-T VariantRecalibrator \\
	-R "${reference}" \\
	--input "${vcf}" \\
	-L "${bed}" \\
	--arg_file TMP/args.list \\
	-mode INDEL \\
	--recal_file "TMP/indel.recal.vcf.gz" \
	--tranches_file "TMP/indel.tranches.txt" \
	-rscriptFile  "TMP/indel.plot.R"


mv -v "TMP/indel.recal.vcf.gz" "${params.prefix}indel.recal.vcf.gz"
mv -v "TMP/indel.recal.vcf.gz.tbi" "${params.prefix}indel.recal.vcf.gz.tbi"
mv -v "TMP/indel.tranches.txt" "${params.prefix}indel.tranches.txt"
"""
}


process APPLY_RECALIBRATION_SNP {
tag "${interval} ${vcf}"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	val(genomeId)
	val(recal)
	val(tranches)
	tuple val(interval),val(vcf)
output:
	tuple val(interval),path("recal.snp.vcf.gz"),emit:output
	path("recal.snp.vcf.gz.tbi")
script:
	def reference = params.genomes[genomeId].fasta

"""
hostname 1>&2
${moduleLoad("java/1.8 bcftools")}
mkdir -p TMP

if ${vcf.endsWith(".bcf")} ; then
	bcftools view --regions "${interval}" --threads ${task.cpus} -o TMP/jeter.vcf.gz -O z "${vcf}"
	bcftools index --threads ${task.cpus}  -t TMP/jeter.vcf.gz
fi

java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
        -T ApplyRecalibration \\
        -R "${reference}" \\
	-L "${interval}" \\
        --input "${vcf.endsWith(".bcf")?"TMP/jeter.vcf.gz":"${vcf}"}" \\
        -mode SNP \\
	--ts_filter_level 99.0 \\
	-recalFile "${recal}" \\
	-tranchesFile "${tranches}" \\
	-o  TMP/recal.snp.vcf.gz

mv -v TMP/recal.* ./

"""
}


process APPLY_RECALIBRATION_INDEL {
tag "${interval} ${vcf}"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	val(genomeId)
	val(recal)
	val(tranches)
	tuple val(interval),val(vcf)
output:
	path("recal.indel.vcf.gz"),emit:output
	path("recal.indel.vcf.gz.tbi")
	
script:
	def reference = params.genomes[genomeId].fasta
"""
hostname 1>&2
${moduleLoad("java/1.8")}
mkdir -p TMP


java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar ${params.gatkjar} \\
        -T ApplyRecalibration \\
        -R "${reference}" \\
	-L "${interval}" \\
        --input "${vcf}" \\
        -mode INDEL \\
	--ts_filter_level 99.0 \\
	-recalFile "${recal}" \\
	-tranchesFile "${tranches}" \\
	-o  TMP/recal.indel.vcf.gz

mv -v TMP/recal.* ./
"""
}

process PLOT_VQSLOD {
afterScript "rm -rf TMP"
tag "${vcf.name}"
input:
	val(meta)
	path(vcf)
output:
	path("${params.prefix}vqsr.pdf"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("r bcftools")} 
mkdir -p TMP

bcftools query -f '%VQSLOD\t%FILTER\t%TYPE\\n' "${vcf}" | sed 's/,OVERLAP\$//' > TMP/jeter.tsv


cat << "__EOF__" > TMP/jeter.R
T1 <- read.table("TMPjeter.tsv",header=FALSE,sep="\t",col.names=c("VQSLOD","FILTER","TYPE"),colClasses=c("numeric","character","character"),stringsAsFactors=FALSE)
head(T1)
pairs = unique(T1[,c("FILTER","TYPE")])
head(pairs)


pdf("${params.prefix}vqsr.pdf")
xminmax = c(-5,5)


for(i in 1:nrow(pairs) ) {
        T2 = T1[T1\$TYPE==pairs[i,"TYPE"] & T1\$FILTER==pairs[i,"FILTER"],]
        d <- density(T2\$VQSLOD)
        plot(d,  xlim = xminmax,
                 ylim = c(0,max(d\$y)),
                xlab = paste(pairs[i,]),
                ylab = "Density",
                main =  paste(pairs[i,]),
                col = ifelse(pairs[i,"FILTER"]=="PASS","green","red")
                ) # plots the results
        }
dev.off()
__EOF__

R --vanilla < TMP/jeter.R

##################################################################################
cat <<- EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot VQSL distribution</entry>
	<entry key="vcf">${vcf}</entry>
</properties>
EOF
"""
}
