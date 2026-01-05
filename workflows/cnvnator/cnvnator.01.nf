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
nextflow.enable.dsl=2

include {isHg19;isHg38;moduleLoad;getVersionCmd;runOnComplete} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include { SCATTER_TO_BED } from '../../subworkflows/picard/picard.scatter2bed.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf'

params.bams=""
params.reference=""
params.help=false
params.ctgRegex="^(chr)?[0-9XY]+\$"
params.conda=""
params.binSizes="5000"
params.prefix=""
params.publishDir=""

workflow {
	cnv_ch = CNVNATOR01(
		params,
		params.reference,
		file(params.bams),
		Channel.from(String.valueOf(params.binSizes).replace("_","").split("[, \t;]+")).filter{S->!S.trim().isEmpty()}
		)
	}

runOnComplete(workflow)


String getGenome(ref) {
	if(isHg38(ref)) return "hg38";
	if(isHg19(ref)) return "hg19";
	log.info("undefined genome "+ref);
	return "";
	}



workflow CNVNATOR01 {
	take:
		meta
		reference
		bams
		bins_ch
	main:
		version_ch = Channel.empty()
		
		one_ctg_ch = ONE_FILE_PER_CONTIG(meta, reference)
		version_ch = version_ch.mix(one_ctg_ch.version)

		all_samples_ch = SAMTOOLS_SAMPLES_01(meta.plus("with_header":false), reference, file("NO_FILE"), bams)
		version_ch= version_ch.mix(all_samples_ch.version)


		each_sample_bam_ch = all_samples_ch.output.
			splitCsv(header:false,sep:'\t').
			map{T->[T[1],T[2]]}
		each_contig  = one_ctg_ch.chroms_list.splitText().map{it.trim()}

		contig_sample_bam_ch = each_contig.combine(each_sample_bam_ch)

		xreads_ch = EXTRACT_READS(meta, reference, contig_sample_bam_ch)
		version_ch = version_ch.mix(xreads_ch.version)

		call_ch = CALL_CNV(meta, reference, one_ctg_ch.chromdir, xreads_ch.output.combine(bins_ch))
		version_ch = version_ch.mix(call_ch.version)

                gaps_ch = SCATTER_TO_BED(["OUTPUT_TYPE":"N","MAX_TO_MERGE":"1"],reference)
		version_ch = version_ch.mix(gaps_ch.version)

		to_zip = Channel.empty()
		merge_ch = MERGE_BIN_SAMPLE(meta, gaps_ch.bed , call_ch.output.map{T->[[T[0],T[2]],[T[1],T[3]] ]}.groupTuple())
		version_ch = version_ch.mix(merge_ch.version)
		to_zip = to_zip.mix(merge_ch.bed.map{T->T[1]})
		to_zip = to_zip.mix(merge_ch.vcf.map{T->T[1]})

		version_ch = MERGE_VERSION(meta,"cnvnator","cnvnator",version_ch.collect())
		to_zip = to_zip.mix(version_ch.version)

		html = VERSION_TO_HTML(params, version_ch)	
		to_zip = to_zip.mix(html.html)


		zip_ch = SIMPLE_ZIP_01([:],to_zip.collect())
	emit:
		version = version_ch
		zip = zip_ch.zip
	}


process ONE_FILE_PER_CONTIG {
tag "${reference}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
output:
	path("CHROMS"),emit:chromdir
	path("CHROMS/chroms.txt"),emit:chroms_list
	path("version.xml"),emit:version
script:
	def ctgRegex = meta.getOrDefault("ctgRegex","^(chr)?[0-9XY]+\$")
"""
hostname 1>&2
${moduleLoad("samtools")}
mkdir TMP

set -x
cut -f 1 "${reference}.fai" | grep -E '${meta.ctgRegex}' | while read C
do
	samtools faidx "${reference}" "\${C}" | tr "\t" " " | cut -f 1 -d ' ' > "TMP/\${C}.fa"
	samtools faidx "TMP/\${C}.fa"
	echo "\${C}" >> TMP/chroms.txt
done

mv TMP CHROMS

##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">generate one file per sequence</entry>
	<entry key="regex">${meta.ctgRegex}</entry>
	<entry key="ref">${reference}</entry>
        <entry key="versions">${getVersionCmd("samtools")}</entry>
</properties>
EOF
"""
}


/* environ 1H15 par contig */
process EXTRACT_READS {
	tag "${sample}/${contig}/${file(bam).name}"
	cache "lenient"
	memory "5g"
	afterScript "rm -rf TMP"
	conda "${meta.conda}/cnvnator"
	when:
		true
	input:
		val(meta)
		val(reference)
		tuple val(contig),val(sample),val(bam)
	output:
		tuple val(contig),val(sample),path("contig.root"),emit:output
		path("version.xml"),emit:version
	script:
		def rgnOpt = "\""+contig+"\""
		def chromOpt = " -chrom "+contig
		def genome = getGenome(reference)
	"""
	hostname 1>&2
	${moduleLoad("samtools")}
	set -o pipefail
	set -x
	mkdir TMP

	if ${bam.endsWith(".cram")} ; then
		samtools view -O BAM -o TMP/jeter.bam -F 3844 --reference "${reference}" "${bam}" "${contig}"
		samtools index TMP/jeter.bam
	fi


	cnvnator \
		-root "TMP/contig.root" \
		${chromOpt} \
		${genome.isEmpty()?"":" -genome "+genome} \
		-tree ${bam.endsWith(".cram")?"TMP/jeter.bam":"'${bam}'"} 1>&2

	mv TMP/contig.root ./

##################
cat << EOF > version.xml
	<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">count reads</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bam">${bam}</entry>
	<entry key="contig">${contig}</entry>
	<entry key="versions">${getVersionCmd("cnvnator samtools")}</entry>
</properties>
EOF
	"""
	}


process CALL_CNV {
	tag "${sample}/${bin}/${contig}"
	memory "5g"
	afterScript "rm -rf TMP"
	conda "${meta.conda}/cnvnator"
	input:
		val(meta)
		val(reference)
		val(contigdir)
		tuple val(contig),val(sample),path(root),val(bin)
	output:
		tuple val(sample),path("cnvnator.bed.gz"),val(bin),path("cnvnator.bcf"),emit:output
		path("version.xml"),emit:version
	script:
		def genome = getGenome(reference)
	"""
	hostname 1>&2
	set -x
	${moduleLoad("bcftools")}
	mkdir TMP
	cp -v "${root}" TMP/tmp.root

	# generate histogram
	cnvnator \
		-root "TMP/tmp.root" \
		-d "${contigdir}" \
		-his ${bin} 1>&2
 
	# calculate statistics
	cnvnator \
		-root "TMP/tmp.root" \
		-d "${contigdir}" \
		-stat ${bin}  1>&2

	# partition
	cnvnator \
		-root "TMP/tmp.root" \
		-d "${contigdir}" \
		-partition ${bin} 1>&2

	# call CNVs
	cnvnator \
		-root "TMP/tmp.root" \
		-d "${contigdir}" \
		-call ${bin} > "TMP/output.tsv" 

	
	cnvnator2VCF.pl \
		-prefix "${bin}" \
		-reference "${genome}" \
		"TMP/output.tsv" \
		"${contigdir}" > TMP/jeter.vcf

	# rename header
	bcftools query -l TMP/jeter.vcf |  awk '{printf("%s\\t${sample}\\n",\$1);}' > TMP/jeter.tsv
	bcftools reheader --fai "${reference}.fai" --samples TMP/jeter.tsv TMP/jeter.vcf > TMP/jeter2.vcf
	mv TMP/jeter2.vcf TMP/jeter.vcf
	
	#sort
	bcftools sort  -o "TMP/jeter.bcf" -O b -T TMP TMP/jeter.vcf
	bcftools index TMP/jeter.bcf


	head TMP/jeter.tsv

	awk -F '\t' '{split(\$2,a,/[:\\-]/);printf("%s\t%d\t%s\t%s\t%s\\n",a[1],int(a[2])-1,a[3],"${sample}",\$0);}' TMP/output.tsv |\
		grep -v 'Number of free parameters' |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		gzip  > TMP/jeter.bed.gz


	mv TMP/jeter.bed.gz ./cnvnator.bed.gz
	mv TMP/jeter.bcf ./cnvnator.bcf
	mv TMP/jeter.bcf.csi ./cnvnator.bcf.csi

##################
cat << EOF > version.xml
	<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">call CNV</entry>
	<entry key="sample">${sample}</entry>
	<entry key="bin">${bin}</entry>
        <entry key="versions">${getVersionCmd("awk bcftools cnvnator")}</entry>
</properties>
EOF
	"""
}


process  MERGE_BIN_SAMPLE {
	tag "${ctgsample[0]}/${ctgsample[1]}/N=${L.size()}"
	conda "${params.conda}/cnvnator"
	when:
		L.size() > 0
	input:
		val(meta)
		path(exclude)
		tuple val(ctgsample),val(L)
	output:
		tuple val(ctgsample),path("${meta.prefix?:""}${ctgsample[0]}.${ctgsample[1]}.cnvnator.bed.gz"),emit:bed
		tuple val(ctgsample),path("${meta.prefix?:""}${ctgsample[0]}.${ctgsample[1]}.cnvnator.bcf"),emit:vcf
		path("version.xml"),emit:version
	script:
		def snx = ctgsample[0]
		def binx = ctgsample[1]
		def prefix= meta.prefix?:""
	"""
	hostname 1>&2
	set -o pipefail
	${moduleLoad("bedtools htslib bcftools")}
	export LC_ALL=C

	gunzip -c  ${L.collect{T->T[0]}.join(" ")} |\
		sort -T . -t '\t'  -k1,1 -k2,2n > jeter.bed


	sort -T . -t '\t'  -k1,1 -k2,2n "${exclude}" > jeter3.bed

	echo "#chrom start end sample CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0" | tr " " "\t" >  jeter2.bed
	bedtools intersect -f 0.95 -r -v -a jeter.bed -b jeter3.bed |\
			sort -T . -t '\t'  -k1,1 -k2,2n >> jeter2.bed
	mv jeter2.bed jeter.bed
	rm jeter3.bed


	bgzip jeter.bed
	mv jeter.bed.gz "${prefix}${snx}.${binx}.cnvnator.bed.gz"
	tabix -p bed "${prefix}${snx}.${binx}.cnvnator.bed.gz"

	bcftools concat --allow-overlaps -O b -o "${prefix}${snx}.${binx}.cnvnator.bcf"  ${L.collect{T->T[1]}.join(" ")} 
	bcftools index "${prefix}${snx}.${binx}.cnvnator.bcf"

##################
cat << EOF > version.xml
	<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">merge cnvnator data for one sample</entry>
	<entry key="sample">${snx}</entry>
	<entry key="bin">${binx}</entry>
        <entry key="versions">${getVersionCmd("bedtools tabix bcftools")}</entry>
</properties>
EOF
	"""
	}

