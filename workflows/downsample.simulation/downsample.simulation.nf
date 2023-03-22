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
nextflow.enable.dsl=2


def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference();

gazoduc.make("bams","NO_FILE").
        description("file containing the path to multiple bam files").
        required().
        existingFile().
        put()

gazoduc.make("mapq",10).
        description("min mapping quality").
        setInteger().
        put()

gazoduc.make("depths","1,2,5,10,15,20").
        description("downsampling depths. comma separated").
        put()


include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {MOSDEPTH_DOWNLOAD_01} from '../../modules/mosdepth/mosdepth.downoad.01.nf'
include {GS_SIMPLE_01} from '../../modules/gs/gs.simple.01.nf'

workflow {
	ch1_ch = DOWNSAMPLE_SIMULATION(params, params.reference, params.bams, params.bed)
	html = VERSION_TO_HTML(params, ch1_ch.version)
}


workflow DOWNSAMPLE_SIMULATION {
take:
	meta
	reference
	bams
	bed
main:
	version_ch  = Channel.empty()

	each_depths= Channel.from(meta.depths.split("[,]")).flatMap().map{T->(T as int)}

	bams_ch = SAMTOOLS_SAMPLES_01(["with_header":false,"allow_multiple_references":false,"allow_duplicate_samples":false], reference, file("NO_FILE"), bams)	
	version_ch = version_ch.mix(bams_ch.version)

	mosdepth_ch = MOSDEPTH_DOWNLOAD_01(meta)
	version_ch = version_ch.mix(mosdepth_ch.version)
	
	sample_bam = bams_ch.output.splitCsv(header:false,sep:'\t').map{T->[T[0],T[2]]}

	cov01_ch = GET_COVERAGE(meta, reference, mosdepth_ch.executable, bed, sample_bam )
	version_ch = version_ch.mix(cov01_ch.version)

	cov02_ch = DOWNSAMPLE(meta, reference, mosdepth_ch.executable, bed, cov01_ch.output.combine(each_depths) )
	version_ch = version_ch.mix(cov02_ch.version)

	gvcf_ch = HAPLOTYPE_CALLER01(meta, reference, bed, 
		sample_bam.map{T->["RAW",T[0],T[1]]}.mix(cov02_ch.output)
		)
	version_ch = version_ch.mix(gvcf_ch.version)


	genotyped_ch = HAPLOTYPE_CALLER02(meta, reference, bed, gvcf_ch.output.map{T->[T[1],T[2]]}.groupTuple())
	version_ch = version_ch.mix(genotyped_ch.version)


	gtc_ch = GENOTYPE_CONCORDANCE(meta, genotyped_ch.output)
	version_ch = version_ch.mix(gtc_ch.version)

	plot_ch=PLOT(meta,gtc_ch.output.collect())
	version_ch = version_ch.mix(plot_ch.version)


	pdf_ch = GS_SIMPLE_01(meta, plot_ch.pdf.mix(gtc_ch.pdf).map{T->["all",T]}.groupTuple())
	version_ch = version_ch.mix(pdf_ch.version)

	version_ch = MERGE_VERSION(meta, "DownSample", "DownSample",version_ch.collect())
emit:
	version = version_ch
	pdf = pdf_ch.output
}

process GET_COVERAGE {
   tag "${sample} / ${file(bam).name}"
   cpus 5
   afterScript "rm *.global.dist.txt *.region.dist.txt *.regions.bed.gz *.regions.bed.gz.csi"
   input:
		val(meta)
		val(reference)
		path(mosdepth)
		path(bed)
		tuple val(sample),val(bam)
   output:
		tuple val(sample),val(bam),path("${sample}.before.mosdepth.summary.txt"),emit:output
		path("version.xml"),emit:version
   script:
   """
   hostname 1>&2

   ${mosdepth.toRealPath()} --no-per-base -t ${task.cpus} --by "${bed}" --fasta "${reference}" \
		--mapq ${meta.mapq} "${sample}.before" "${bam}"
	
###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">get bam depth</entry>
        <entry key="sample">${sample}</entry>
        <entry key="versions"></entry>
</properties>
EOF
   """
}

process DOWNSAMPLE {
   tag "${sample} / ${file(bam).name} DP=${depth}"
   memory '5g'
   afterScript "rm -rf TMP"
   input:
	val(meta)
	val(reference)
	path(mosdepth)
	path(bed)
	tuple val(sample),val(bam),path(summary),val(depth)
   output:
	tuple val(depth),val(sample),path("${sample}.DP${depth}.bam"),emit:output
	path("version.xml"),emit:version
   script:
	def dps = String.format("%03d",depth)
   """
   hostname 1>&2
   ${moduleLoad("samtools gatk/0.0.0")}
   mkdir -p TMP

   cut -f 1  "${bed}" | uniq | sort | uniq | while read C
   do
	FRAC=`awk -F '\t' -vC=\$C '(\$1==sprintf("%s_region",C)) {D=(\$4*1.0);if(D==0) D=1E-6; F=(${depth}*1.0)/D ; if(F>0 && F<1.0) print F;}' '${summary}'  `

	echo "${sample} \${C} FRAC=\${FRAC}" >&2 

	awk -F '\t' -vC=\$C '(\$1==C)' "${bed}" > TMP/jeter.bed


	if [ "\${FRAC}" != "" ] ; then

   		samtools view -M -L TMP/jeter.bed --threads ${task.cpus} -q ${meta.mapq} -F 3844 -u --reference ${reference} ${bam}  |\
			samtools view -s \${FRAC} -O BAM -o "TMP/_chunck.\${C}.bam"

	else

   		samtools view -M -L TMP/jeter.bed --threads ${task.cpus} -q ${meta.mapq} -F 3844 -u --reference ${reference} \
			-O BAM -o "TMP/_chunck.\${C}.bam" ${bam} "\${C}" 
      
	fi
   done

   samtools merge --threads ${task.cpus} --reference ${params.reference} TMP/merged.bam TMP/_chunck.*.bam

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" AddOrReplaceReadGroups \
	     I=TMP/merged.bam \
	     O="${sample}.DP${depth}.bam" \
	     RGID=${sample}.DP${dps} \
	     RGLB=${sample}.DP${dps} \
	     RGSM=${sample}.DP${dps} \
	     RGPL=ILLUMINA \
	     RGPU=unit1

   samtools index -@ ${task.cpus} "${sample}.DP${depth}.bam"

   ${mosdepth.toRealPath()} --no-per-base -t ${task.cpus} --by ${bed} --fasta ${reference} \
		--mapq ${params.mapq} TMP/${sample}.after.DP${depth} ${sample}.DP${depth}.bam

   mv "TMP/${sample}.after.DP${depth}.mosdepth.summary.txt" ./

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">downsample BAM</entry>
        <entry key="sample">${sample}</entry>
        <entry key="depth">${depth}</entry>
        <entry key="versions">${getVersionCmd("samtools")}</entry>
</properties>
EOF
   """
   }


process HAPLOTYPE_CALLER01 {
   tag "${sample} ${bam.name}"
   memory '5g'
   afterScript "rm -rf TMP"
   input:
	val(meta)
        val(reference)
	path(bed)
	tuple val(depth),val(sample),path(bam)
   output:
	tuple val(depth),val(sample),path("genotype.g.vcf.gz"),emit:output
	path("version.xml"),emit:version
   script:
   """ 
   hostname
   set -o pipefail
   ${moduleLoad("gatk/0.0.0")}
   mkdir -p TMP


   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
	     -I "${bam.toRealPath()}" \
	     -ERC GVCF \
	     -L "${bed}" \
	     -R "${reference}" \
	     --minimum-mapping-quality  ${meta.mapq} \
	     --verbosity ERROR \
	     -O "genotype.g.vcf.gz"
	

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">downsample BAM</entry>
        <entry key="sample">${sample}</entry>
        <entry key="bam">${bam}</entry>
        <entry key="versions">${getVersionCmd("gatk")}</entry>
</properties>
EOF

   """
}



process HAPLOTYPE_CALLER02 {
   tag "${sample} N=${L.size()}"
   cache "lenient"
   memory '5g'
   afterScript "rm -rf TMP"
   input:
	val(meta)
        val(reference)
	path(bed)
	tuple val(sample),val(L)
   output:
	tuple val(sample),path("${sample}.vcf.gz"),emit:output
	path("version.xml"),emit:version
   script:
   """ 
   hostname
   set -o pipefail
   ${moduleLoad("gatk/0.0.0 bcftools")}
   mkdir -p TMP

cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" CombineGVCFs \
	-R "${reference}" \
	-L "${bed}" \
	--verbosity ERROR \
	-V TMP/jeter.list \
	-O TMP/combined.g.vcf.gz

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeGVCFs \
        -R "${reference}" \
        -L "${bed}" \
        --verbosity ERROR \
        -V TMP/combined.g.vcf.gz \
        -O TMP/jeter.vcf.gz \
	-G StandardAnnotation -G AS_StandardAnnotation

   bcftools view -O b -o "${sample}.vcf.gz" TMP/jeter.vcf.gz
   bcftools index -t "${sample}.vcf.gz"
	

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">call BAMs</entry>
        <entry key="sample">${sample}</entry>
        <entry key="bams">${L.size()}</entry>
        <entry key="versions">${getVersionCmd("gatk")}</entry>
</properties>
EOF
   """
}

process GENOTYPE_CONCORDANCE {
   tag "${sample} N=${vcf.name}"
   cache "lenient"
   memory '5g'
   //afterScript "rm -rf TMP"
   input:
        val(meta)
        tuple val(sample),path(vcf)
   output:
	path("concordances.txt"),emit:output
	path("${sample}.GQ.pdf"),emit:pdf
        path("version.xml"),emit:version
   script:
   """
   hostname 1>&2
   set -o pipefail
   ${moduleLoad("gatk/0.0.0 bcftools r")}
   mkdir -p TMP

   bcftools query -l "${vcf}" | awk '(\$1!="${sample}")' | while read S
   do

   gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeConcordance \
	--TRUTH_VCF "${vcf.toRealPath()}" \
	--TRUTH_SAMPLE "${sample}" \
	--CALL_VCF "${vcf.toRealPath()}" \
	--CALL_SAMPLE "\${S}" \
	--O "\${S}"

    grep -E '^(SNP|INDEL)' "\${S}.genotype_concordance_summary_metrics"  |\
	grep -E '^(SNP|INDEL)' | cut -f 1,13 |\
	awk -vS=\$S '{printf("%s\t%s\\n", gensub(/.*\\.DP[0]*/,"\\\\1","g",S) ,\$0);}' >> concordances.txt
   done


   bcftools query -l "${vcf}" | while read S
   do
      bcftools view -O u --samples "\${S}" "${vcf}" |\
	bcftools query -f '[%GQ\\n]' | awk '\$1!="."' | sort -n >> "TMP/\${S}.dist"
   done


echo -n "filenames <-c(" > TMP/jeter.R

find TMP -type f -name "*.dist" | sort | awk '{printf("\\"%s\\"\\n",\$0);}' |paste -s -d, >> TMP/jeter.R

cat << '__EOF__' >> TMP/jeter.R
)
n <- 1
colors <- rainbow(length(filenames))
pdf("${sample}.GQ.pdf")
list1 <- list()
list2 <- c()
for(filename in filenames) {
	X <- scan(filename, comment.char = ".")
	list2 <- append(list2,gsub("${sample}.DP","",gsub(".dist","",basename(filename))))
	if(n==1) {
		plot(density(X),type='l',main="${sample} Genotype Quality",xlab="GQ",ylab="Fraction",col=colors[n])
	} else	{		
		lines(density(X),col=colors[n])
		}
	n <- n+1
}

legend("top", legend=list2,text.col=colors)
dev.off()
__EOF__

R --vanilla < TMP/jeter.R

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Genotype Concordance</entry>
        <entry key="sample">${sample}</entry>
        <entry key="vcf">${vcf}</entry>
        <entry key="versions">${getVersionCmd("gatk bcftools")}</entry>
</properties>
EOF
"""
}

process PLOT {
tag "${files.size()}"
input:
	val(meta)
	val(files)
output:
	path("${meta.prefix?:""}concordance.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("r")}


cat ${files.join(" ")} > jeter.tsv

cat << '__EOF__' > jeter.R

T1 <- read.table("jeter.tsv",sep="\t",header=FALSE)
depths <- sort(as.integer(unique(T1\$V1)), decreasing = TRUE)
types <- unique(T1\$V2)

list1 <- list()
list2 <- c()
list3 <- c()

for(type in types) {
for(depth in depths) {

	T2 <-T1[as.integer(T1\$V1)==depth & T1\$V2==type,]\$V3
	list1[[length(list1)+1]] <- as.double(T2)
	list2 <- append(list2,depth)
	list3 <- append(list3,ifelse(type=="SNP","blue","yellow"))
}
}


pdf("${meta.prefix?:""}concordance.pdf")

boxplot(list1,names=list2,main="${meta.prefix?:""}concordance",xlab="depth",ylab="type",col=list3,las=2)

legend("topright",legend=c("SNP","INDEL"),fill=c("blue","yellow"))

dev.off()

__EOF__

R --vanilla < jeter.R

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">plot Genotype Concordance</entry>
        <entry key="n.data">${files.size()}</entry>
        <entry key="versions">${getVersionCmd("R")}</entry>
</properties>
EOF
"""
}

