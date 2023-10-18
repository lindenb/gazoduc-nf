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
include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {SCATTER_TO_BED} from '../../subworkflows/picard/picard.scatter2bed.nf'
include {MOSDEPTH_DOWNLOAD_01} from '../../modules/mosdepth/mosdepth.downoad.01.nf'
include {MOSDEPTH_RUN_01} from '../../modules/mosdepth/mosdepth.run.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {moduleLoad;getVersionCmd} from '../../modules/utils/functions.nf'


workflow MOSDEPTH_BAMS_01 {
	take:
		meta
		reference
		bams
		bed
	main:
		version_ch  = Channel.empty()
		bams_ch = SAMTOOLS_SAMPLES_01(meta.plus(["with_header":true,"allow_multiple_references":false,"allow_duplicate_samples":false]),reference,file("NO_FILE"),bams)
		version_ch = version_ch.mix(bams_ch.version)

		if(bed.name.equals("NO_FILE")) {
			acgt_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"ACGT","MAX_TO_MERGE":"1"],reference)
                        version_ch = version_ch.mix(acgt_ch.version)
                        bed2 = acgt_ch.bed
			}
		else
			{
			bed2 = Channel.fromPath(bed) // rechanged ? changed 20230517 ... for graphtyper ?? Channel.fromPath(bed)
			}

		ch1 = bams_ch.output.splitCsv(header:true,sep:'\t').combine(bed2).map{T->T[0].plus([
			"mapq": (meta.mapq?:1),
			"bed": T[1]
			])}
		mosdepth_ch = MOSDEPTH_DOWNLOAD_01(meta)
		version_ch = version_ch.mix(mosdepth_ch.version)	
	
		ch2 = MOSDEPTH_RUN_01(meta, mosdepth_ch.executable,ch1)
		version_ch = version_ch.mix(ch2.version)

		merge_ch = MERGE_MOSDEPTH_SUMMARY(meta,ch2.summary.map{T->T[0].sample+"\t"+T[0].bam+"\t"+T[0].reference+"\t"+T[1]}.collect())
		version_ch = version_ch.mix(merge_ch.version)

		file_list_ch = COLLECT_TO_FILE_01([:],ch2.regiondist.map{T->T[1]}.collect())
		version_ch = version_ch.mix(file_list_ch.version)

		regdist_ch = ZIP_REGION_DIST(meta, file_list_ch.output)
		version_ch = version_ch.mix(regdist_ch.version)

		multiqc_ch = MULTIQC_01(meta,ch2.regiondist.map{T->T[1]}.collect())
		version_ch = version_ch.mix(multiqc_ch.version)

		pdf_ch = Channel.empty()
		plot_ch = PLOT_IT(meta, merge_ch.output)
		version_ch = version_ch.mix(plot_ch.version)
		pdf_ch = pdf_ch.mix(plot_ch.coverage)
		pdf_ch = pdf_ch.mix(plot_ch.coverage_region)


		ch3 = MERGE_VERSION(meta, "Mosdepth", "Mosdepth", version_ch.collect())
	emit:
		version = ch3
		summary = merge_ch.output
		globaldist = ch2.globaldist
		perbase = ch2.perbase
		pdf  = pdf_ch.collect()
		region_dist_zip = regdist_ch.zip
		multiqc = multiqc_ch.zip
	}

process MERGE_MOSDEPTH_SUMMARY {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(L)
output:
	path("${meta.prefix?:""}summary.tsv"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("r")}
set -o pipefail

mkdir TMP

cat << EOF > TMP/jeter1.txt
${L.join("\n")}
EOF

echo "sample\tbam\treference\tcoverage\tcoverage_region" > TMP/jeter2.txt

cat TMP/jeter1.txt | while read SN BAM REF SUMMARY
do
	echo -n "\${SN}\t\${BAM}\t\${REF}" >> TMP/jeter2.txt
	awk -F '\t' '(\$1=="total") {printf("\t%s",\$4);}' "\${SUMMARY}" >> TMP/jeter2.txt
	awk -F '\t' 'BEGIN{C="";} (\$1=="total") {if(C=="") {C=\$4;}} (\$1=="total_region") {C=\$4;} END{printf("\t%s",C);}' "\${SUMMARY}" >> TMP/jeter2.txt
	echo >> TMP/jeter2.txt
done

mv TMP/jeter2.txt "${meta.prefix?:""}summary.tsv"

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge mosdepth summary</entry>
	<entry key="count">${L.size()}</entry>
</properties>
EOF
"""
}

process ZIP_REGION_DIST {
tag "${files.name}"
input:
	val(meta)
	path(files)
output:
	path("${meta.prefix?:""}regions.dist.zip"),emit:zip
	path("version.xml"),emit:version
script:
"""
hostname 1>&2

zip -9@j "${meta.prefix?:""}regions.dist.zip" < "${files}"

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge mosdepth region_dist</entry>
        <entry key="input">${files}</entry>
</properties>
EOF
"""
}

process PLOT_IT {
tag "${summary.name}"
input:
	val(meta)
	path(summary)
output:
	path("${meta.prefix?:""}coverage.pdf"),emit:coverage
	path("${meta.prefix?:""}coverage_region.pdf"),optional:true,emit:coverage_region
	path("version.xml"),emit:version
script:
	def prefix = meta.prefix?:""
"""
hostname 1>&2
${moduleLoad("r")}

cat << EOF > jeter.R
T1 <- read.table(file="${summary}",sep="\\t",header=TRUE)
T2 <- as.numeric(T1\\\$coverage)

pdf("${prefix}coverage.pdf")
boxplot(T2 ,ylim=c(0,max(T2)+1), main="${prefix}mosdepth",xlab="Sample",ylab="coverage")
dev.off()


T1 <- T1[T1\\\$coverage_region!="",]

if (nrow(T1)>0) {
T2 <- as.numeric(T1\\\$coverage_region)

pdf("${prefix}coverage_region.pdf")
boxplot(T2 ,ylim=c(0,max(T2)+1), main="${prefix}mosdepth in BED",xlab="Sample",ylab="coverage")
dev.off()
}


EOF

R --no-save < jeter.R

rm jeter.R

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot mosdepth summary</entry>
	<entry key="versions">${getVersionCmd("r")}</entry>
</properties>
EOF
"""
}
