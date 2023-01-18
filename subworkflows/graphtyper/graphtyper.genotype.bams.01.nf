/*

Copyright (c) 2023 Pierre Lindenbaum

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
include {isBlank;moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'
include {MOSDEPTH_BAMS_01} from '../../subworkflows/mosdepth/mosdepth.01.nf'
include {COLLECT_TO_FILE_01 as COLLECT_TO_FILE_X1; COLLECT_TO_FILE_01 as COLLECT_TO_FILE_X2} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../../subworkflows/bcftools/bcftools.concat.01.nf'
include {SAMTOOLS_DEPTH_01} from '../../subworkflows/samtools/samtools.depth.01.nf'
include {SAMTOOLS_SAMPLES_01} from '../samtools/samtools.samples.01.nf'



workflow GRAPHTYPER_GENOTYPE_BAMS_01 {
	take:
		meta
		reference
		bams
		beds
		sample2depth
	main:
		version_ch = Channel.empty()
		
		
		concat_bed_ch = CONCAT_BEDS(meta, beds)
		version_ch = version_ch.mix(concat_bed_ch.version)
		
		
		/** file SAMPLE <-> DEPTH provided */
		if(!sample2depth.name.equals("NO_FILE")) {
			bams_ch = SAMTOOLS_SAMPLES_01(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false],reference,file("NO_FILE"),bams)
			version_ch = version_ch.mix(bams_ch.version)

			join1_ch = JOIN_WITH_DEPTH(meta,bams_ch.output,sample2depth)

			summary = join1_ch.output
			}
		else if((meta.depth_method?:"mosdepth").equals("samtoolsdepth")) {
			stdp_ch = SAMTOOLS_DEPTH_01(meta,reference,file("NO_FILE"),bams, concat_bed_ch.bed)
			version_ch = version_ch.mix(stdp_ch.version)
			
			convert_ch = CONVERT_TO_MOSDEPTH_SUMMARY(meta,stdp_ch.output)

			summary= convert_ch.output
			}
		else 
			{
			mosdepth_ch = MOSDEPTH_BAMS_01(meta, reference, bams, concat_bed_ch.bed)
			version_ch = version_ch.mix(mosdepth_ch.version)
			summary= mosdepth_ch.summary
			}
		
		x1_ch = COVERAGE_DIVIDE_READLENGTH(meta,summary)
		version_ch = version_ch.mix(x1_ch.version)
		

		executable_ch = GRAPHTYPER_DOWNLOAD_01(meta)
		version_ch = version_ch.mix(executable_ch.version)

		
		each_bed_ch = Channel.fromPath(beds).splitText().map{it.trim()}
	
		x2_ch = GRAPHTYPER_GENOTYPE_01(meta, executable_ch.executable,x1_ch.output.combine(each_bed_ch).map{T->[
			"bams":T[0],
			"avg_cov_by_readlen":T[1],
			"bed":T[2],
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

process CONCAT_BEDS {
executor "local"
input:
	val(meta)
	path(beds)
output:
	path("concat.bed"),emit:bed
	path("version.xml"),emit:version
script:
"""
set -o pipefail
${moduleLoad("bedtools")}
xargs -a "${beds}" cut -f1,2,3 |\
	sort -T . -t '\t' -k1,1 -k2,2n |\
	bedtools merge > concat.bed

test -s concat.bed

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">concatenate bed files</entry>
	<entry key="version">${getVersionCmd("bedtools")}</entry>
</properties>
EOF
"""
}

process JOIN_WITH_DEPTH {
executor "local"
input:
      	val(meta)
        path(samplebam)
        path(sample2depth)
output:
        path("summary.tsv"),emit:output
script:
"""
echo -e "sample\tbam\tref\tcov\tcovr" > jeter.tsv

join -t '\t' -1 1 -2 1 -o '1.1,1.3,1.4,2.2,2.2' \
	<(sort -t '\t' -k1,1 "${samplebam}" ) \
	<(sort -t '\t' -k1,1 "${sample2depth}" )  >> jeter.tsv


mv jeter.tsv summary.tsv
"""
}


process CONVERT_TO_MOSDEPTH_SUMMARY {
executor "local"
input:
	val(meta)
	path(depths)
output:
	path("summary.tsv"),emit:output
script:
"""
echo -e "sample\tbam\tref\tcov\tcovr" > jeter.tsv
awk -F '\t' '/^#/ {next;} {printf("%s\t%s\t%s\t%s\t%s\\n",\$1,\$2,\$3,\$6,\$6);}' '${depths}' >> jeter.tsv


mv jeter.tsv summary.tsv
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
