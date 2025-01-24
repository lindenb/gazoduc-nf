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

include {moduleLoad;runOnComplete} from '../../../modules/utils/functions.nf'
//include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
//include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
//include {JVARKIT_VCF_TO_BED_01} from '../jvarkit/jvarkit.vcf2bed.01.nf'


workflow  {

		reference = Channel.of(file(params.fasta),file(params.fai),file(params.dict)).collect()

		if(params.vcf.endsWith(".list")) {
			vcf_ch = Channel.fromPath(params.vcf).splitText().map{it.trim()}
			}
		else
			{
			vcf_ch = Channel.of(params.vcf)
			}
		vcf_ch = vcf_ch.map{[it,it+(it.endsWith(".bcf")?".csi":".tbi")]}
		
		
		vcf2bed_ch = VCF_TO_BED(vcf_ch)

		vcf2inter_ch = VCF_TO_INTERVALS(vcf2bed_ch.output.splitCsv(header:false,sep:'\t').map{it.flatten()}, file(params.bed))

		//vcf2inter_ch.splitCsv(header:false,sep:'\t').view()
		snpeff_ch = DOWNLOAD_SNPEFF_DB()


		annot_ch = ANNOT(
			reference,
			snpeff_ch.output ,
			file(params.bed),
			vcf2inter_ch.output.splitCsv(header:false,sep:'\t').map{it.flatten()}.filter{!it[0].equals(".")}
			)


		CONCAT(annot_ch.output.flatten().collect())
		
	}


process VCF_TO_BED {
label "process_quick"
input:
	tuple path(vcf),path(vcfidx)
output:
	tuple path("vcf.bed"),path(vcf),path(vcfidx),emit:output
script:
"""
module load bcftools bedtools
set -o pipefail

bcftools index -s "${vcf}" |\
	awk -F '\t' '{printf("%s\t0\t%s\\n",\$1,\$2);}' |\
	sort -t '\t' -k1,1 -k2,2n -T . > vcf.bed
"""
}


process VCF_TO_INTERVALS {
tag "${contig} ${start0} ${end} ${vcf} ${vcfidx} ${bed}"
label "process_quick"
afterScript "rm -rf TMP"
input:
	tuple val(contig),val(start0),val(end),path(vcf),path(vcfidx)
	path(bed)
output:
	tuple path("intervals.bed"),path(vcf),path(vcfidx),emit:output
script:
"""
mkdir -p TMP
module load bedtools bcftools jvarkit

echo -e "${contig}\t${start0}\t${end}" > TMP/jeter1.bed

if ${!bed.name.equals("NO_FILE")}
then
	sort -t '\t' -k1,1 -k2,2n -T TMP ${bed} > TMP/jeter2.bed
	bedtools intersect -a TMP/jeter1.bed -b TMP/jeter2.bed > TMP/jeter3.bed
	mv TMP/jeter3.bed TMP/jeter1.bed
fi

if test  -s TMP/jeter1.bed
then

	bcftools view -G --regions-file TMP/jeter1.bed  "${vcf}" |\\
		java -jar -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP \${JVARKIT_DIST}/jvarkit.jar vcf2intervals --distance ${params.vcf2interval_distance} --bed |\\
		cut -f1,2,3 > TMP/intervals.bed

	mv TMP/intervals.bed ./
else
	echo -e ".\t0\t1" > intervals.bed
fi

"""
}



process  DOWNLOAD_SNPEFF_DB {
tag "${params.snpeff_db}"
label "process_short"
output:
        path("SNPEFF"),emit:output
script:
"""
module load snpEff/5.2.c
mkdir -p SNPEFFX TMP
java -Xmx${task.memory.giga}g  -Djava.io.tmpdir=TMP -jar \${SNPEFF_JAR} \\
	download \\
	-dataDir "\${PWD}/SNPEFFX"  \\
	'${params.snpeff_db}'
test -s SNPEFFX/*/snpEffectPredictor.bin
mv SNPEFFX SNPEFF
"""
}


process ANNOT {
tag "${contig}:${start0}-${end}"
label "process_quick"
afterScript "rm -rf TMP"
memory "5g"
input:
	path(reference)
	path(snpeffdir)
	path(bed)
	tuple val(contig),val(start0),val(end),path(vcf),path(vcfidx)
output:
	path("${contig}_${start0}_${end}.*"),emit:output
script:
	def gnomadAF = params.gnomadAF
	def gnomadPop = params.gnomadPop
	def soacn = params.soacn
	def soft_filters = params.soft_filters.toLowerCase().split("[, ]");
	def soft_filter_gnomad = soft_filters.contains("gnomad")
	def soft_filter_so = soft_filters.contains("so")
"""
hostname 1>&2
module load bedtools bcftools jvarkit snpEff/5.2.c
mkdir -p TMP


echo -e "${contig}\t${start0}\t${end}" > TMP/jeter.bed

if ${!bed.name.equals("NO_FILE")}
then
	awk -F '\t' '(\$1=="${contig}")' "${bed}" | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/jeter2.bed
	bedtools intersect -a TMP/jeter.bed -b TMP/jeter2.bed > TMP/jeter3.bed
	mv TMP/jeter3.bed  TMP/jeter.bed

fi

bcftools view --regions-file TMP/jeter.bed -O u  -o TMP/jeter1.bcf "${vcf}"


bcftools view TMP/jeter1.bcf |\\
java -jar -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP \${SNPEFF_JAR} eff \\
	-dataDir "\${PWD}/${snpeffdir}" \\
	-nodownload -noNextProt -noMotif -noInteraction -noLog -noStats -chr chr -i vcf -o vcf "${params.snpeff_db}" > TMP/jeter1.vcf


if ${!soacn.isEmpty()} ; then
	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcffilterso \\
		${soft_filter_so?"--filterout BAD_SO":""} \\
		--acn "${soacn}" TMP/jeter1.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter1.vcf
fi

java -Xmx${task.memory.giga}G  -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/jvarkit.jar vcfgnomad  --bufferSize 10000 \
		--gnomad "${params.gnomad}" \
		--fields "${gnomadPop}" \
		--max-af "${gnomadAF}" TMP/jeter1.vcf > TMP/jeter2.vcf

mv TMP/jeter2.vcf TMP/jeter1.vcf

if ${!soft_filter_gnomad}
then
	bcftools view -e 'FILTER ~ "GNOMAD_GENOME_InbreedingCoeff" || FILTER ~ "GNOMAD_GENOME_AC0" || FILTER  ~ "GNOMAD_GENOME_AS_VQSR" || FILTER ~ "GNOMAD_GENOME_BAD_AF"' TMP/jeter1.vcf > TMP/jeter2.vcf

	mv TMP/jeter2.vcf TMP/jeter1.vcf
fi


if ${params.containsKey("cadd") && !params.cadd.isEmpty()}
then

	java -Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP  -jar \${JVARKIT_DIST}/jvarkit.jar vcfcadd \
                        --tabix "${params.cadd}" TMP/jeter1.vcf > TMP/jeter2.vcf
        mv TMP/jeter2.vcf TMP/jeter1.vcf
fi

bcftools view TMP/jeter1.vcf -O b9 -o "${contig}_${start0}_${end}.annot.bcf"

bcftools index -f "${contig}_${start0}_${end}.annot.bcf"
"""
}



process CONCAT {
label "process_short"
input:
	path("VCF/*")
output:
	path("annot.*"),emit:output
script:
"""
module load bcftools
find VCF  \\( -name "*.vcf.gz" -o -name "*.bcf" \\) > jeter.list
test -s jeter.list
bcftools concat -a -d all --file-list jeter.list -O b9 -o annot.bcf
rm jeter.list
bcftools index -f annot.bcf
"""
}
