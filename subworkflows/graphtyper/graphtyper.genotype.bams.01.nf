/*

Copyright (c) 2022 Pierre Lindenbaum

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
include {GRAPHTYPER_DOWNLOAD_01} from '../../modules/graphtyper/graphtyper.download.01.nf'
include {GRAPHTYPER_GENOTYPE_01} from '../../modules/graphtyper/graphtyper.genotype.01.nf'
include {SCATTER_TO_BED} from '../../subworkflows/picard/picard.scatter2bed.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {isBlank;moduleLoad} from '../../modules/utils/functions.nf'
include {MOSDEPTH_BAMS_01} from '../../subworkflows/mosdepth/mosdepth.01.nf'
include {COLLECT_TO_FILE_01 as COLLECT_TO_FILE_X1; COLLECT_TO_FILE_01 as COLLECT_TO_FILE_X2} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'
include {SAMTOOLS_SAMPLES01} from '../../modules/samtools/samtools.samples.01.nf'

workflow GRAPHTYPER_GENOTYPE_BAMS_01 {
	take:
		meta
		reference
		bams
		bed
	main:
		version_ch = Channel.empty()
		
		each_interval_ch = bed.splitCsv(header:false,sep:'\t').
			map{T->T[0]+":"+((T[1] as int)+1)+"-"+T[2]}

		if((meta.depth_method?:"mosdepth").equals("samtoolsdepth")) {
			snbam_ch = SAMTOOLS_SAMPLES01(meta, referenc, bams)
			version_ch = version_ch.mix(snbal_ch.version)

			each_sn_bam_ch = snbam_ch.output.splitCsv(header:false,sep:'\t')

			stdepth_ch = SAMTOOLS_DEPTH(meta, reference, each_sn_bam_ch)
			version_ch = version_ch.mix(stdepth_ch.version)

			COLLECT_TO_FILE_X2(
			
			}
		else	
			{
			gaps_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"1"],reference)
			version_ch = version_ch.mix(gaps_ch.version)		

			mosdepth_ch = MOSDEPTH_BAMS_01(meta, reference, bams, gaps_ch.bed)
			version_ch = version_ch.mix(mosdepth_ch.version)
			}
		
		x1_ch = COVERAGE_DIVIDE_READLENGTH(meta,mosdepth_ch.summary)
		version_ch = version_ch.mix(x1_ch.version)


		executable_ch = GRAPHTYPER_DOWNLOAD_01(meta)
		version_ch = version_ch.mix(executable_ch.version)

		x2_ch = GRAPHTYPER_GENOTYPE_01(meta, executable_ch.executable,x1_ch.output.combine(each_interval_ch).map{T->[
			"bams":T[0],
			"avg_cov_by_readlen":T[1],
			"interval":T[2],
			"reference":reference
			]})
		version_ch = version_ch.mix(x2_ch.version)

		x3_ch = COLLECT_TO_FILE_X1([:],x2_ch.output.map{T->T[1]}.collect())
		version_ch = version_ch.mix(x3_ch.version)


		x4_ch = BCFTOOLS_CONCAT_01(meta,x3_ch.output)
		version_ch = version_ch.mix(x4_ch.version)

		version_ch = MERGE_VERSION(meta, "graptyper", "genotype bams with graphtyper", version_ch.collect())
	emit:
		version = version_ch
		vcf = x4_ch.vcf
	}

process SAMTOOLS_DEPTH {
input:
	val(meta)
	val(reference)
	path(bed)
	tuple val(sample),val(bam)
script:
	def mapq=meta.mapq?:30
"""
hostname 1>&2
${moduleLoad("samtools")}
set -o pipefail

samtools view -F 3844 -T "${reference}" --min-MQ "${mapq}" -M -L "${bed}" -u -O BAM "${bam}" |\
samtools depth -a -Q "${mapq}" |\
awk -F '\t' 'BEGIN{N=0;T=0.0;printf("#sample\tbam\t\cov\tcov\n");}	{T+=int(\$3);N++;} END {printf("${sample}\t${bam}\t${reference}\t%f\t%f\\n",T/N,T/N);}' > ${sample}.depth.tsv

##############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">use samtools depth to get depth in bed</entry>
        <entry key="bed">${bed}</entry>
</properties>
EOF
"""
}

process COVERAGE_DIVIDE_READLENGTH {
tag "${summary.name}"
input:
	val(meta)
	path(summary)
output:
	tuple path("bams.txt"),path("avg_cov_by_readlen.txt"),emit:output
	path("version.xml"),emit:version
script:
	def n_reads = meta.n_reads?:1000
"""
${moduleLoad("samtools")}
tail -n+2 "${summary}" | while read SN BAM REF COV COVR
do
	samtools view -F 3844 -T "\${REF}" "\${BAM}" | head -n "${n_reads}" |\
		awk -F '\t' -vCOVR=\${COVR} 'BEGIN{T=0.0;N=0;} {N++;T+=length(\$10)} END{print COVR/(N==0?100:T/N);}' >> avg_cov_by_readlen.txt
	echo "\${BAM}" >> bams.txt
done

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">create list of bams and list of coverage/read-length</entry>
</properties>
EOF
"""
}
