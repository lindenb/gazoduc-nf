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
nextflow.enable.dsl=2


/*
inspiration:
	* https://genome.cshlp.org/content/32/3/569.full
	* https://api.firecloud.org/ga4gh/v1/tools/mitochondria:AlignAndCall/versions/23/plain-WDL/descriptor
	* https://github.com/grbot/varcall/blob/5e5f744ca891415b6997c1072d4fcf3743bf92d7/mt-calling/main.nf
	* https://api.firecloud.org/ga4gh/v1/tools/mitochondria:AlignmentPipeline/versions/1/plain-WDL/descriptor

*/

/** path to indexed fasta reference */
params.reference = ""
params.references = "NO_FILE"
params.prefix = ""
params.publishDir = ""
params.bams = "NO_FILE"
params.mt_accession = "NC_012920.1"

include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {getVersionCmd;moduleLoad;runOnComplete} from '../../modules/utils/functions.nf'
include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'

def helpMessage() {
  log.info"""
## About

call MT genome

## Author

${params.rsrc.author}

## Options

  * --reference (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) path to multiple bams files [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--vcf in.vcf.gz \\
	--bed /path/to/in.bed
```

## Workflow

![workflow](./workflow.svg)
  
## See also


"""
}


if( params.help ) {
    helpMessage()
    exit 0
}


workflow {
	ch1 = MT_CALLING(params, params.reference, file(params.references), file(params.bams))
	html = VERSION_TO_HTML(params,ch1.version)
	//PUBLISH(ch1.version,html.html,ch1.summary,ch1.pdf.collect())
	}

workflow MT_CALLING {
	take:
		meta
		reference
		references
		bams
	main:
		version_ch = Channel.empty()

		mt_ch = MT_REFERENCE(meta)
		version_ch = version_ch.mix(mt_ch.version)

		shifted_ch = SHIFTED_REFERENCE(meta, mt_ch.reference)
		version_ch = version_ch.mix(shifted_ch.version)

		bams_ch = SAMTOOLS_SAMPLES_01(meta.plus(["with_header":true,"allow_multiple_references":true,"allow_duplicate_samples":true]),reference,references,bams)
		version_ch = version_ch.mix(bams_ch.version)


		each_bam = bams_ch.output.splitCsv(header:true,sep:'\t')

		sn_ch = EXTRACT_MT_BAM(meta, mt_ch.reference, shifted_ch.reference, each_bam)
		version_ch = version_ch.mix(sn_ch.version)
		
		merge_ch = MERGE_BCFS(meta,sn_ch.vcf.collect())
		version_ch = version_ch.mix(merge_ch.version)

		pdf_ch = MERGE_PDFS(meta,sn_ch.pdf.collect())
		version_ch = version_ch.mix(pdf_ch.version)

		version_ch = MERGE_VERSION(meta, "Mitochondrial Calling", "Mitrochondrial Call", version_ch.collect())
	emit:
		version= version_ch
		vcf = merge_ch.vcf
		index = merge_ch.index
		pdf = pdf_ch.pdf
	}

process MT_REFERENCE {
input:
	val(meta)
output:
	path("chrM.fa"),emit:reference
	path("version.xml"),emit:version
script:
	def acn = meta.mt_accession
"""
hostname 1>&2
${moduleLoad("samtools bwa")}
wget -O chrM.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acn}&rettype=fasta&retmode=text"
samtools faidx chrM.fa
samtools dict -A -o chrM.dict chrM.fa
bwa index -a is chrM.fa

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">download and index mitochondrial genome</entry>
	<entry key="acn">${acn}</entry>
	<entry key="versions">${getVersionCmd("bwa samtools wget")}</entry>
</properties>
EOF
"""
}


process SHIFTED_REFERENCE {
input:
	val(meta)
	path(reference)
output:
	path("chrM_shifted.fa"),emit:reference
	path("version.xml"),emit:version
script:
	def ref = reference.toRealPath()
"""
hostname 1>&2
${moduleLoad("samtools bwa")}
set -o pipefail

# find half length to shift reference
MT=`awk -F '\t' '(NR==1){printf("%s\\n",\$1);}' "${ref}.fai"`
LN=`awk -F '\t' '(NR==1){printf("%s\\n",\$2);}' "${ref}.fai"`
SHIFT=\$((\${LN}/2))
test ! -z "\${MT}"

head -n 1 "${ref}" > tmp.fa
samtools faidx "${ref}" "\${MT}:\$((\${SHIFT}+1))-\${LN}" | grep -v "^>" | tr -d '\\n' >> tmp.fa
samtools faidx "${ref}" "\${MT}:1-\${SHIFT}" | grep -v "^>" | tr -d '\\n' >> tmp.fa

fold -w 100 < tmp.fa > chrM_shifted.fa
rm tmp.fa

samtools faidx chrM_shifted.fa
samtools dict -A -o chrM_shifted.dict chrM_shifted.fa
bwa index -a is chrM_shifted.fa

###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">shift mitochondrial genome</entry>
	<entry key="reference">${ref}</entry>
	<entry key="shit">\${SHIFT}</entry>
	<entry key="versions">${getVersionCmd("bwa samtools awk")}</entry>
</properties>
EOF
"""
}


process EXTRACT_MT_BAM {
tag "${row.new_sample} ${file(row.bam).name}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(meta)
	path(reference)
	path(shifted_reference)
	val(row)
output:
	path("${row.sample}.genotyped.bcf"),emit:vcf
	path("${row.sample}.coverage.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
	def ref1  = reference.toRealPath()
	def ref2  = shifted_reference.toRealPath()
"""
hostname 1>&2
${moduleLoad("samtools bwa gatk/0.0.0 bcftools R")}
set -o pipefail
set -x
mkdir TMP

cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">remap and call bam</entry>
	<entry key="reference.out">${ref1}</entry>
	<entry key="shifted.reference.out">${ref2}</entry>
	<entry key="sample">${row.new_sample}</entry>
	<entry key="bam">${row.bam}</entry>
	<entry key="reference.in">${row.reference}</entry>
	<entry key="versions">${getVersionCmd("samtools bwa gatk awk bcftools")}</entry>
EOF

MT=`cut -f 1 "${row.reference}.fai" | awk '\$1 ~ /^(chr)?M(T)?\$/' | head -n 1`
echo "#MT is (\${MT})" 1>&2
test ! -z "\${MT}"

## extract mitochondrial reads
samtools view -F '3852' --expr 'rname == rnext' -T "${row.reference}" --uncompressed -O BAM "${row.bam}" "\${MT}" |\
samtools collate -f -O -u --no-PG --reference "${row.reference}" - TMP/tmp.collate |\
samtools fastq -N -1 TMP/jeter.R1.fq.gz -2 TMP/jeter.R2.fq.gz -s TMP/jeter.Rx.fq.gz -0 TMP/jeter.R0.fq.gz -n

cat << EOF >> version.xml
	 <entry key="fastq.read.count">\$(gunzip -c TMP/jeter.R1.fq.gz | paste - - - - | wc -l)</entry>
EOF


# name and size mitochondrion
MT=`awk -F '\t' '(NR==1){printf("%s\\n",\$1);}' "${ref1}.fai"`
LN=`awk -F '\t' '(NR==1){printf("%s\\n",\$2);}' "${ref1}.fai"`
SHIFT=\$((\${LN}/2))


i=1

# map on reference and shifted reference
for REF in "${ref1}" "${ref2}"
do

	# only call in that area, exclude boundaries
	awk -F '\t' '(NR==1){L=int(\$2);X=500;printf("%s\t%d\t%d\\n",\$1,X,L-X);}' "\${REF}.fai" > TMP/jeter.bed
	

	## Y use clipping for sup align
	bwa mem -Y -R '@RG\\tID:${row.new_sample}\\tSM:${row.new_sample}\\tLB:${row.new_sample}\\tCN:NantesBird\\tPL:ILLUMINA' "\${REF}" TMP/jeter.R1.fq.gz TMP/jeter.R2.fq.gz |\
	samtools view -F '3844' --expr '(rname == rnext) || (!flag.unmap && flag.munmap)' -o TMP/jeter.bam -O BAM -


	cat <<- EOF >> version.xml
		<entry key="remapped.read.count.\${i}">\$(samtools view -c TMP/jeter.bam)</entry>
	EOF


	# remove duplicates
	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \
	      MarkDuplicates \
	      INPUT=TMP/jeter.bam \
	      OUTPUT=TMP/jeter2.bam \
	      METRICS_FILE=TMP/markdup.metrics \
	      VALIDATION_STRINGENCY=SILENT \
	      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
	      ASSUME_SORT_ORDER=queryname \
	      REMOVE_SEQUENCING_DUPLICATES=true \
	      CLEAR_DT="false" \
	      ADD_PG_TAG_TO_READS=false

	mv TMP/jeter2.bam TMP/jeter.bam

	# sort BAM
	samtools sort -o TMP/jeter2.bam -O BAM -T TMP/tmpx TMP/jeter.bam
	mv TMP/jeter2.bam TMP/jeter.bam
	samtools  index TMP/jeter.bam

	# save depth  for plot
	samtools depth -a -r "\${MT}" TMP/jeter.bam | cut -f 2,3 > TMP/jeter.depth.tsv

	if test "\$i"  -eq 1 ; then
		mv TMP/jeter.depth.tsv "TMP/jeter.depth.\${i}.tsv"
	else
		awk -F '\t' -vSHIFT=\${SHIFT} -vLN=\${LN} '{OFS="\t";P=int(\$1);P=P-(SHIFT+1);if(P<1) P=LN+P;\$1=P;print;}' TMP/jeter.depth.tsv |\
			sort -t '\t' -T TMP -k1,1n > "TMP/jeter.depth.\${i}.tsv"
		rm TMP/jeter.depth.tsv 

		## merge both depth
		paste TMP/jeter.depth.1.tsv TMP/jeter.depth.2.tsv |\
			awk -F '\t' 'BEGIN{printf("POS\tT1\tT2\\n");} (\$1==\$3) {printf("%s\t%s\t%s\\n",\$1,\$2,\$4);}' > TMP/jeter.depth.tsv
		rm TMP/jeter.depth.1.tsv TMP/jeter.depth.2.tsv

		# plot coverage
		cat <<- 'EOF' | R --vanilla
		T1<-read.table("TMP/jeter.depth.tsv",header=T,sep="\t",stringsAsFactor=F,colClasses=c("integer","integer","integer"))
		
		runmedk<-10
		
		pdf("${row.sample}.coverage.pdf",width=20)
		plot(runmed(T1\$T1,runmedk),type="l",
			main="Coverage on MT genome for ${row.sample}",
			sub="${row.bam}",
			col="blue",
			ylab="depth",
			xlab="position",
			xlim=c(1,max(T1\$POS)),
			ylim=c(0,max(1,max(T1\$T1,T1\$T2)))
			)
		
		par(new=TRUE)
		
		plot(runmed(T1\$T2,runmedk),type="l",
			main="",
			sub="",
			xlab="",
			ylab="",
			col="green",
			xlim=c(1,max(T1\$POS)),
			ylim=c(0,max(1,max(T1\$T1,T1\$T2)))
			)
		legend("topleft",legend=c("REF","SHITED REF"),
			fill=c("blue","green"))
		dev.off();
		EOF

	fi


	# collect metrics
	gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" \
	      CollectWgsMetrics \
	      INPUT=TMP/jeter.bam \
	      VALIDATION_STRINGENCY=SILENT \
	      REFERENCE_SEQUENCE=\${REF} \
	      OUTPUT=wgs.metrics.\${i}.txt \
	      USE_FAST_ALGORITHM=true \
	      INCLUDE_BQ_HISTOGRAM=true \
	      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.\${i}.txt

	# call mutect
	gatk  --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" Mutect2 \
		-R "\${REF}" \
		-L TMP/jeter.bed \
		--mitochondria-mode \
		--annotation StrandBiasBySample \
		--max-reads-per-alignment-start  75 \
		--max-mnp-distance 0 \
		-I TMP/jeter.bam \
		-O TMP/jeter.\${i}.vcf.gz


	if test "\$i"  -eq 1 ; then
		bcftools view -O v  TMP/jeter.1.vcf.gz |\
			awk '/^#CHROM/ {print("##INFO=<ID=PPOS,Number=1,Type=Integer,Description=\\"POS in ${ref2} BEFORE SHIFT\\">\\n");} {print;}' |\
			bcftools view -O b -o TMP/jeter.1.bcf
			
	else
		# use the header of the original MT genome
		bcftools view --header-only TMP/jeter.1.bcf > TMP/jeter2.vcf
		# shift all POS to the original MT genome
		bcftools view --no-header TMP/jeter.2.vcf.gz |\
			awk -F '\t' -vSHIFT=\${SHIFT} -vLN=\${LN} '{OFS="\t";\$8=sprintf("%s;PPOS=%s",\$8,\$2);P=int(\$2) - (SHIFT+1); if(P<1) {P=LN+P;} ; \$2=P;print;}' >> TMP/jeter2.vcf
		bcftools sort -T TMP -O b -o TMP/jeter.2.bcf TMP/jeter2.vcf
		rm TMP/jeter2.vcf
	fi

	bcftools index TMP/jeter.\${i}.bcf

i=\$((i+1))
done



# merge VCF mapped on ref and shifted ref
bcftools concat  --allow-overlaps --remove-duplicates -O z -o TMP/jeter.vcf.gz TMP/jeter.1.bcf TMP/jeter.2.bcf
rm TMP/jeter.1.bcf TMP/jeter.2.bcf
bcftools index --tbi --force TMP/jeter.vcf.gz

# check no problem with shift
bcftools norm --check-ref e --fasta-ref "${ref1}" TMP/jeter.vcf.gz > /dev/null


# merge mutect stats
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" MergeMutectStats \
        --stats TMP/jeter.1.vcf.gz.stats --stats TMP/jeter.2.vcf.gz.stats \
        -O TMP/merged.stats

# cleanup Mutect2
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" FilterMutectCalls \
	-V TMP/jeter.vcf.gz \
        -R "${ref1}" \
        -O TMP/jeter2.vcf.gz \
        --stats "TMP/merged.stats" \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0 

mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
mv TMP/jeter2.vcf.gz.tbi TMP/jeter.vcf.gz.tbi

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" VariantFiltration \
	-V TMP/jeter.vcf.gz \
        -O TMP/jeter2.vcf.gz \
        --apply-allele-specific-filters 

mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
mv TMP/jeter2.vcf.gz.tbi TMP/jeter.vcf.gz.tbi

# bug mutect2 https://github.com/broadinstitute/gatk/issues/6857
bcftools view TMP/jeter.vcf.gz |\
	sed '/ID=AS_FilterStatus/s/Number=A/Number=./'  |\
	bcftools view -O b -o TMP/jeter.bcf -


mv TMP/jeter.bcf "${row.sample}.genotyped.bcf"
bcftools index "${row.sample}.genotyped.bcf"

cat << EOF >> version.xml
</properties>
EOF
"""
}



process MERGE_BCFS {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
memory "3g"
input:
        val(meta)
        val(L)
output:
        path("${meta.prefix?:""}genotyped.bcf"),emit:vcf
        path("${meta.prefix?:""}genotyped.bcf.csi"),emit:index
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("bcftools")}
set -o pipefail
set -x

cat << EOF > jeter.list
${L.join("\n")}
EOF

bcftools merge --missing-to-ref --file-list jeter.list -O b -o "${meta.prefix?:""}genotyped.bcf"
bcftools index "${meta.prefix?:""}genotyped.bcf"

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge bcfs</entry>
        <entry key="count">${L.size()}</entry>
        <entry key="versions">${getVersionCmd("bcftools")}</entry>
</properties>
EOF
"""
}

process MERGE_PDFS {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
memory "3g"
input:
        val(meta)
        val(L)
output:
        path("${meta.prefix?:""}coverage.pdf"),emit:pdf
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -o pipefail
set -x

cat << EOF | awk -F '/' '{printf("%s\t%s\\n",\$NF,\$0);}' | sort -t '\t' -T . -k1,1 | cut -f2 > jeter.list
${L.join("\n")}
EOF

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="${meta.prefix?:""}coverage.pdf" @jeter.list


cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">merge pdf</entry>
        <entry key="count">${L.size()}</entry>
        <entry key="versions">${getVersionCmd("gs")}</entry>
</properties>
EOF
"""
}



process PUBLISH {
publishDir "${params.publishDir}" , mode: 'copy', overwrite: true
input:
	path(version)
	path(html)
	path(summary)
	path(pdfs)
output:
	path("${params.prefix}mosdepth.zip")
when:
	!params.getOrDefault("publishDir","").trim().isEmpty()
script:
"""
zip -j "${params.prefix}mosdepth.zip" ${version} ${html} ${summary} ${pdfs.join(" ")}
"""
}

runOnComplete(workflow);

