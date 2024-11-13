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
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run my_pipeline --input input_file.csv")
   exit 0
}
// validate parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


include {BCFTOOLS_CONCAT} from '../../subworkflows/bcftools/bcftools.concat.02.nf'
include {moduleLoad } from '../../modules/utils/functions.nf'



workflow {
		if(params.vcf.endsWith(".list")) {
			in_vcf_ch = Channel.fromPath(params.vcf).
				splitText().
				map{it.trim()}.
				map{file(it)}
			}
		else
			{
			in_vcf_ch = Channel.fromPath(params.vcf)
			}
		in_vcf_ch = in_vcf_ch.map{[it, file(it.name.endsWith(".bcf")?""+it+".csi":""+it+".tbi")]}.collect()
		vcf2bed_ch = BCF_TO_VCF(in_vcf_ch)

		compile_ch= COMPILE_MINIKIT()

		each_contig_ch = vcf2bed_ch.contigs.splitText().map{it.trim()}

		vcf2rgn_ch = VCF_TO_INTERVALS(vcf2bed_ch.output.combine(each_contig_ch))



		genome_ch = [file(params.fasta) , file(params.fasta+".fai"), file(""+file(params.fasta).getParent()+"/"+file(params.fasta).getBaseName()+".dict" ) ]	



		recal_snp_ch = RECALIBRATE_SNP( genome_ch, vcf2bed_ch.output , compile_ch.output, file(params.recal_snps) )

		recal_indel_ch = RECALIBRATE_INDEL(genome_ch, vcf2bed_ch.output ,compile_ch.output, file(params.recal_indels))


		rgn_vcf_ch = vcf2rgn_ch.output.splitCsv(sep:'\t',header:false).
			map{ it[0]+":"+((it[1] as int)+1)+"-"+it[2] }.
			combine(recal_snp_ch.output).
			combine(recal_indel_ch.output)


		apply_snp_ch = APPLY_RECALIBRATION_SNP(genome_ch, in_vcf_ch, rgn_vcf_ch)
		apply_indel_ch = APPLY_RECALIBRATION_INDEL(genome_ch , apply_snp_ch.output)

		if(params.gather_by.equals("chrom") || params.gather_by.equals("chromosome") || params.gather_by.equals("contig")) {
			chx = each_contig_ch.combine(apply_indel_ch.flatten().collect().toList())
			chx.view()
			//GATHER_BY_CONTIG(chx)
			}
		else
			{
			GATHER_VCFS(apply_indel_ch.flatten().collect())
			}
		
	}


process BCF_TO_VCF {
label "process_short"
afterScript 'rm -rf  TMP'
input:
	path("VCFS/*")
output:
	tuple path("reformat.vcf.gz"),path("reformat.vcf.gz.tbi"),emit:output
	path("chroms.txt"),emit:contigs
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
mkdir -p TMP
find VCFS -name "*.vcf.gz" -o -name "*.bcf" > TMP/jeter.list


bcftools concat --threads ${task.cpus} --drop-genotypes -a -O z -o TMP/reformat.vcf.gz --file-list TMP/jeter.list
bcftools index  --threads ${task.cpus}  --tbi TMP/reformat.vcf.gz 

mv TMP/reformat.* ./
	
bcftools index --stats reformat.vcf.gz | cut -f 1 > chroms.txt

"""
}


process COMPILE_MINIKIT {
label "process_short"
afterScript "rm -rf TMP"
output:
	path("minikit.jar"),emit:output
script:
"""
${moduleLoad("java/17.0.10")}

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

process VCF_TO_INTERVALS {
tag "${contig} ${vcf}"
label "process_short"
input:
	tuple path(vcf),path(idx),val(contig)
output:
	path("intervals.bed"),emit:output
script:
	def distance="10Mb"
	def min_distance=100
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit")}

	bcftools view -G "${vcf}" "${contig}" |\
	java -Xmx${task.memory.giga}G -jar \${JVARKIT_DIST}/jvarkit.jar vcf2intervals \
		--bed \
		--distance "${distance}" \
		--min-distance "${min_distance}" > intervals.bed

"""
}

process RECALIBRATE_SNP {
tag "${vcf.name}"
label "process_medium"
afterScript 'rm -rf  TMP'
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(vcf),path(vcfidx)
	path(minikit)
	path(recal_snps)
output:
	tuple path("snp.recal.vcf.gz"),path("snp.recal.vcf.gz.tbi"),path("snp.tranches.txt"),emit:output
	path("snps.plot.R"),emit:R
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}
set -x
mkdir -p TMP

test -s "${recal_snps}"

bcftools view --type snps -G "${vcf}" | java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum > TMP/args.list
test -s TMP/args.list


cat ${recal_snps} >> TMP/args.list


gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -XX:-UsePerfData" VariantRecalibrator  \\
	-R "${fasta}" \\
	-V "${vcf}" \
	--arguments_file TMP/args.list \\
	-mode SNP \\
	-O "snp.recal.vcf.gz" \\
	--tranches-file "snp.tranches.txt" \\
        --dont-run-rscript \\
	--rscript-file  "snps.plot.R"

"""
}


process RECALIBRATE_INDEL {
tag "${vcf}"
label "process_medium"
cpus 5
memory "15g"
afterScript 'rm -rf  TMP'
input:
	tuple path(fasta),path(fai),path(dict)
	tuple path(vcf),path(vcfidx)
	path(minikit)
	path(recal_indels)
output:
	tuple path("indel.recal.vcf.gz"),path("indel.recal.vcf.gz.tbi"),path("indel.tranches.txt"),emit:output
	path("indels.plot.R"),emit:R
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}
mkdir -p TMP
set -x

test -s "${recal_indels}"

bcftools view --type indels -G "${vcf}" |\\
	java -jar ${minikit} --an ExcessHet,FS,InbreedingCoeff,QD,ReadPosRankSum,SOR,MQ,MQRankSum > TMP/args.list

test -s TMP/args.list


cat "${recal_indels}" >> TMP/args.list

# 20200630 add AS https://gatk.broadinstitute.org/hc/en-us/articles/360036510892-VariantRecalibrator
# finalement AS ne semble pas marcher avec indel

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantRecalibrator  \\
	-R "${fasta}" \\
	-V "${vcf}" \\
	--arguments_file TMP/args.list \\
	-mode INDEL \\
	-O "indel.recal.vcf.gz" \\
	--tranches-file "indel.tranches.txt" \\
        --dont-run-rscript \\
	--rscript-file indels.plot.R
"""
}

process APPLY_RECALIBRATION_SNP {
tag "${interval}"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	tuple path(fasta),path(fai),path(dict)
	path("VCFS/*")
	tuple val(interval),
		path(recal_snp_vcf),path(recal_snp_vcf_idx),path(recal_snp_tranches),
		path(recal_indel_vcf),path(recal_indel_vcf_idx),path(recal_indel_tranches)
output:
	tuple val(interval),
		path("recal.snp.vcf.gz"),path("recal.snp.vcf.gz.tbi"),
		path(recal_indel_vcf),path(recal_indel_vcf_idx),path(recal_indel_tranches),emit:output
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r bcftools")}
mkdir -p TMP

find VCFS/ -name "*.vcf.gz" -o -name "*.bcf" > TMP/vcfs.list
bcftools concat -a --regions "${interval}" --file-list TMP/vcfs.list -O z -o TMP/jeter.vcf.gz
bcftools index -t -f TMP/jeter.vcf.gz

gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" ApplyVQSR  \
        -R "${fasta}" \
        -V "TMP/jeter.vcf.gz"  \
        -mode SNP \
        -L "${interval}" \
        --truth-sensitivity-filter-level 99.0 \
        --recal-file "${recal_snp_vcf}" \
        --tranches-file "${recal_snp_tranches}" \
        -O "TMP/jeter2.vcf.gz"

mv TMP/jeter2.vcf.gz recal.snp.vcf.gz
mv TMP/jeter2.vcf.gz.tbi recal.snp.vcf.gz.tbi

"""
}


process APPLY_RECALIBRATION_INDEL {
tag "${interval}"
label "process_medium"
afterScript 'rm -rf  TMP'
memory '10g'
input:
	tuple path(fasta),path(fai),path(dict)
	tuple val(interval),path(vcf),path(vcf_idx),
		path(recal_indel_vcf),path(recal_indel_vcf_idx),path(recal_indel_tranches)
output:
	tuple path("*.recal.vcf.gz"),path("*.recal.vcf.gz.tbi"),emit:output
script:
	def suffix = interval.replaceAll("[^A-Za-z_0-9]+","_");
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 r")}
mkdir -p TMP

gatk  --java-options "-Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP" ApplyVQSR  \
	-R "${fasta}" \
	-V "${vcf}"  \
	-mode INDEL \
	-L "${interval}" \
	--truth-sensitivity-filter-level 99.0 \\
	--recal-file "${recal_indel_vcf}" \\
	--tranches-file "${recal_indel_tranches}" \\
	-O "TMP/jeter.vcf.gz"

mv TMP/jeter.vcf.gz "${suffix}.recal.vcf.gz"
mv TMP/jeter.vcf.gz.tbi "${suffix}.recal.vcf.gz.tbi"
"""
}



process GATHER_VCFS {
label "process_medium"
afterScript "rm -rf TMP"
cpus 10
memory "10G"
input:
	path("VCFS/*")
output:
	tuple path("vqsr.bcf"),path("vqsr.bcf.csi"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 bcftools")} 
mkdir -p TMP

find VCFS/ -name "*.vcf.gz" > TMP/jeter.list

gatk  --java-options "-Xmx${task.memory.giga}g -XX:-UsePerfData -Djava.io.tmpdir=TMP" GatherVcfs  \\
	--INPUT TMP/jeter.list \\
	--REORDER_INPUT_BY_FIRST_VARIANT \\
	--OUTPUT TMP/jeter.vcf.gz	

bcftools view --threads "${task.cpus}" -O b9 -o TMP/jeter.bcf TMP/jeter.vcf.gz
bcftools index --force --threads "${task.cpus}" TMP/jeter.bcf

mv TMP/jeter.bcf vqsr.bcf
mv TMP/jeter.bcf.csi vqsr.bcf.csi
"""
}

process GATHER_BY_CONTIG {
label "process_medium"
tag "${contig}"
afterScript "rm -rf TMP"
cpus 10
memory "10G"
input:
	tuple val(contig),path("VCFS/*")
output:
	tuple path("${contig}.vqsr.bcf"),path("${contig}.vqsr.bcf.csi"),emit:output
script:
"""
hostname 1>&2
${moduleLoad("bcftools")} 
mkdir -p TMP

find VCFS/ -name "*.vcf.gz" | awk -F '/' '{printf("%s,%s\\n",\$NF,\$0);}' | LC_ALL=C sort -t, -k1,1V -T TMP | cut -d, -f2  > TMP/jeter.list
test -s TMP/jeter.list

bcftools concat -a --regions "${contig}" --threads "${task.cpus}" --file-list TMP/jeter.list -O b9 -o TMP/jeter.bcf
bcftools index --force --threads "${task.cpus}" TMP/jeter.bcf

mv TMP/jeter.bcf ${contig}.vqsr.bcf
mv TMP/jeter.bcf.csi ${contig}.vqsr.bcf.csi
"""
}
