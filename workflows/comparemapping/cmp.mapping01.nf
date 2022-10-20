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

/** path to indexed fasta reference */
params.reference1 = ""
params.reference2 = ""
params.bams = "NO_FILE"
/** vcf on ref 1 */
params.vcf1 = "NO_FILE"
params.help = false
params.publishDir = ""
params.prefix = ""

include {SAMTOOLS_SAMPLES_01} from '../../subworkflows/samtools/samtools.samples.01.nf'
include {VERSION_TO_HTML} from '../../modules/version/version2html.nf'
include {moduleLoad;runOnComplete;parseBoolean;getVersionCmd} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {GET_LIFTOVER_CHAINS as LIFT1; GET_LIFTOVER_CHAINS as LIFT2} from '../../modules/ucsc/liftover.chains.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'

def helpMessage() {
  log.info"""
## About

apply mosdepth to a set of bams.

## Author

${params.rsrc.author}

## Options

  * --reference1 (fasta) ${params.rsrc.reference} [REQUIRED]
  * --reference2 (fasta) ${params.rsrc.reference} [REQUIRED]
  * --bams (file) one file containing the paths to the BAM/CRAM [REQUIRED]
  * --publishDir (dir) Save output in this directory
  * --prefix (string) files prefix. default: ""

## Usage

```
nextflow -C ../../confs/cluster.cfg  run -resume workflow1.nf \\
	--publishDir output \\
	--prefix "analysis." \\
	--reference /path/to/reference.fasta \\
	--bams /path/to/bams.list \\
	--bed /path/to/file.bed
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
	ch1 = COMPARE_MAPPING(params,params.reference1,params.reference2,file(params.bams),params.vcf1)
	//PUBLISH(ch1.zip)
	}

runOnComplete(workflow)

workflow COMPARE_MAPPING {
	take:
		meta
		reference1
		reference2
		bams
		vcf1
	main:
		version_ch  = Channel.empty()

		references = REFERENCES(meta,reference1,reference2)

		bams_ch = SAMTOOLS_SAMPLES_01(["with_header":false,"allow_multiple_references":true,"allow_duplicate_samples":true],reference1,references.output,bams)
		version_ch = version_ch.mix(bams_ch.version)
	
		chain12_ch = LIFT2([:],[name:"ref1ToRef2","reference1":reference1,"reference2":reference2])
		version_ch = version_ch.mix(chain12_ch.version)

		chain21_ch = LIFT1([:],[name:"ref2ToRef1","reference1":reference2,"reference2":reference1])
		version_ch = version_ch.mix(chain21_ch.version)


		vcf_ch = GET_VCF1(meta, vcf1, reference1,reference2,chain12_ch.chain,chain21_ch.chain)
		version_ch = version_ch.mix(vcf_ch.version)

		join_ch=JOIN(meta, reference1, reference2, bams_ch.output, vcf_ch.variants_ref1, vcf_ch.variants_ref2,chain21_ch.chain)
		version_ch = version_ch.mix(join_ch.version)

		call_ch= CALL_VARIANT(meta,join_ch.output.splitCsv(header:true,sep:'\t'))
		version_ch = version_ch.mix(call_ch.version)


		join2_ch=JOIN2(meta,reference1, reference2, call_ch.output.map{T->T.join("\t")}.collect())
		version_ch = version_ch.mix(join2_ch.version)

		gt_ch = GT_CONCORDANCE(meta,join2_ch.output.splitCsv(header:true,sep:'\t'))
		version_ch = version_ch.mix(gt_ch.version)


		col_ch = COLLECT_TO_FILE_01([:],gt_ch.output.collect())
		version_ch = version_ch.mix(col_ch.version)

		plot_ch = PLOT(meta,col_ch.output)
		version_ch = version_ch.mix(plot_ch.version)

		version_ch = MERGE_VERSION(meta, "CompareMapping", "Compare Mapping",version_ch.collect())
		

		html = VERSION_TO_HTML(meta,version_ch.version)
		
		zip_ch = ZIP(meta,col_ch.output, version_ch.version, html.html, plot_ch.pdf)
	emit:
		zip = zip_ch.zip
		pdf = plot_ch.pdf
		version= version_ch
	}

process REFERENCES {
executor "local"
input:
	val(meta)
	val(reference1)
	val(reference2)
output:
	path("references.txt"),emit:output
script:
"""
echo "${reference1}" >  references.txt
echo "${reference2}" >> references.txt
"""
}

process  GET_VCF1 {
tag "${vcf1}"
afterScript "rm -rf TMP"
memory "10g"
input:
	val(meta)
	val(vcf1)
	val(reference1)
	val(reference2)
	path(ref1_to_ref2_chain)
	path(ref2_to_ref1_chain)
output:
	path("variants.ref1.vcf.gz"),emit:variants_ref1
	path("variants.ref2.vcf.gz"),emit:variants_ref2
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("jvarkit bcftools gatk/0.0.0")}
set -x

mkdir -p TMP
wget -q -O - "${vcf1}" |\
	bcftools annotate --set-id '+%CHROM-%POS-%REF-%FIRST_ALT' -O v |\
	java -Djava.io.tmpdir=TMP -jar \${JVARKIT_DIST}/vcfsetdict.jar -n SKIP -R "${reference1}"  |\
	awk -F '\t' '\$1 ~ /^#/ || \$1 ~ /^(chr)?[0-9]+\$/' |\
	bcftools sort -T TMP -O z -o TMP/jeter1.vcf.gz

bcftools index -t TMP/jeter1.vcf.gz

bcftools query -f '%ID\\n' TMP/jeter1.vcf.gz | sort -T TMP  | uniq -u > TMP/id1.txt
test -s TMP/id1.txt

## LIFTOVER TO REF2
gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" LiftoverVcf \
	--CHAIN ${ref1_to_ref2_chain} \
	--INPUT TMP/jeter1.vcf.gz \
	--OUTPUT TMP/jeter2.vcf.gz \
	--REFERENCE_SEQUENCE ${reference2} \
	--REJECT TMP/reject.vcf.gz \
	--LOG_FAILED_INTERVALS false \
	--WRITE_ORIGINAL_ALLELES false \
	--WRITE_ORIGINAL_POSITION false


bcftools query -f '%ID\\n' TMP/jeter2.vcf.gz | sort -T TMP | uniq -u > TMP/id2.txt

comm -12 TMP/id1.txt TMP/id2.txt > TMP/id12.txt
test -s TMP/id12.txt

bcftools view -i 'ID=@TMP/id12.txt' -O u TMP/jeter1.vcf.gz  |\
bcftools sort  --max-mem ${task.memory.giga}G -T TMP -O z -o "TMP/variants.ref1.vcf.gz"
bcftools index -t TMP/variants.ref1.vcf.gz

bcftools view -i 'ID=@TMP/id12.txt' -O u TMP/jeter2.vcf.gz  |\
bcftools sort  --max-mem ${task.memory.giga}G -T TMP -O z -o "TMP/variants.ref2.vcf.gz"
bcftools index -t TMP/variants.ref2.vcf.gz

mv TMP/variants.ref1.vcf.gz ./
mv TMP/variants.ref1.vcf.gz.tbi ./
mv TMP/variants.ref2.vcf.gz ./
mv TMP/variants.ref2.vcf.gz.tbi ./


###############################################################################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">get VCF on reference1</entry>
        <entry key="vcf1">${vcf1}</entry>
        <entry key="versions">${getVersionCmd("wget bcftools gatk jvarkit/vcfsetdict")}</entry>
</properties>
EOF
"""
}

process JOIN {
executor "local"
afterScript "rm -f sn1.txt sn2.txt sn12.txt jeter1.txt jeter2.txt jeter3.txt"
input:
	val(meta)
	val(reference1)
	val(reference2)
	path(sn_bam)
	path(vcf1)
	path(vcf2)
	path(chain21)
output:
	path("output.tsv"),emit:output	
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -x

awk -F '\t' '(\$4=="${reference1}") {printf("%s\t${vcf1.toRealPath()}\tNO_FILE\tNO_FILE\\n",\$0);}' '${sn_bam}' | sort -T . -k1,1 -t '\t' > jeter1.txt
cut -f1 jeter1.txt | sort | uniq -u > sn1.txt
test -s jeter1.txt

awk -F '\t' '(\$4=="${reference2}") {printf("%s\t${vcf2.toRealPath()}\t${chain21.toRealPath()}\t${reference1}\\n",\$0);}' '${sn_bam}' | sort -T . -k1,1 -t '\t' > jeter2.txt
cut -f1 jeter2.txt | sort | uniq -u > sn2.txt
test -s jeter2.txt

comm -12 sn1.txt sn2.txt > sn12.txt

join -t '\t' -1 1 -2 1 sn12.txt jeter1.txt > jeter3.txt && mv jeter3.txt jeter1.txt && test -s jeter1.txt
join -t '\t' -1 1 -2 1 sn12.txt jeter2.txt > jeter3.txt && mv jeter3.txt jeter2.txt && test -s jeter2.txt

echo -e "sample\tnew_sample\tbam\treference\tvcf\tchain\tdest_reference" > output.tsv
cat jeter2.txt jeter1.txt >> output.tsv


hostname 1>&2
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">Dispatch calling per sample</entry>
</properties>
EOF
"""
}


process CALL_VARIANT {
tag "${row.sample} ${file(row.bam).name} ${file(row.reference).name}"
memory "5g"
cpus 2
afterScript "rm -rf TMP"
input:
	val(meta)
	val(row)
output:
	tuple val("${row.sample}"),val("${row.reference}"),path("genotyped.vcf.gz"),emit:output
	path("version.xml"),emit:version
script:
	def mapq = meta.mapq?:"30"
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0 bcftools")}

mkdir TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" HaplotypeCaller \
	-I "${row.bam}" \
	-L "${row.vcf}" \
	-G StandardHCAnnotation -G StandardAnnotation \
	--do-not-run-physical-phasing \
	--alleles "${row.vcf}"  \
	--force-call-filtered-alleles \
	--minimum-mapping-quality ${mapq} \
	--dbsnp "${row.vcf}" \
	-R "${row.reference}" \
	-O "TMP/jeter.vcf.gz"

# be sure there is an ID
bcftools view --known -O z -o TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz

bcftools index --force -t TMP/jeter.vcf.gz

if ${!row.chain.equals("NO_FILE")} ; then

## change name for future join

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP"  RenameSampleInVcf \
        --INPUT TMP/jeter.vcf.gz \
        --OUTPUT TMP/jeter2.vcf.gz \
	--NEW_SAMPLE_NAME  "${row.sample}_2"

mv TMP/jeter2.vcf.gz TMP/jeter.vcf.gz
bcftools index --force -t TMP/jeter.vcf.gz


## LIFTOVER TO REF2

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" LiftoverVcf \
	--CHAIN "${row.chain}" \
	--INPUT TMP/jeter.vcf.gz \
	--OUTPUT TMP/jeter2.vcf.gz \
	--REFERENCE_SEQUENCE "${row.dest_reference}" \
	--REJECT TMP/reject.vcf.gz \
	--LOG_FAILED_INTERVALS false \
	--WRITE_ORIGINAL_ALLELES false \
	--WRITE_ORIGINAL_POSITION false

	mv TMP/jeter2.vcf.gz genotyped.vcf.gz
else
	mv TMP/jeter.vcf.gz genotyped.vcf.gz

fi

bcftools index -t genotyped.vcf.gz

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">call sample and liftover to ref1 if needed</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="bam">${row.bam}</entry>
        <entry key="ref">${row.reference}</entry>
        <entry key="vcf">${row.vcf}</entry>
        <entry key="versions">${getVersionCmd("bcftools gatk")}</entry>
</properties>
EOF
"""
}


process JOIN2 {
tag "N=${L.size()}"
executor "local"
afterScript "rm -f jeter.tsv jeter1.tsv jeter2.tsv"
input:
      	val(meta)
        val(reference1)
        val(reference2)
        val(L)
output:
       	path("output.tsv"),emit:output
        path("version.xml"),emit:version
script:
"""
hostname 1>&2
set -x

cat << EOF > jeter.tsv
${L.join("\n")}
EOF

awk -F '\t' '(\$2=="${reference1}")' jeter.tsv | cut -f1,3 | sort -T . -k1,1 -t '\t' > jeter1.tsv
test -s jeter1.tsv
awk -F '\t' '(\$2=="${reference2}")' jeter.tsv | cut -f1,3 | sort -T . -k1,1 -t '\t' > jeter2.tsv
test -s jeter2.tsv

echo -e "sample\tvcf1\tvcf2" > output.tsv
join -t '\t' -1 1 -2 1 -o '1.1,1.2,2.2' jeter1.tsv jeter2.tsv >> output.tsv

cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join calling on REF1 and REF2 with sample</entry>
        <entry key="versions">${getVersionCmd("awk")}</entry>
</properties>
EOF
"""
}

process GT_CONCORDANCE {
tag "${row.sample}"
afterScript "rm -rf TMP"
memory "3g"
input:
	val(meta)
	val(row)
output:
	path("${meta.prefix?:""}${row.sample}.genotype_concordance_summary_metrics"),emit:output
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("gatk/0.0.0")}

mkdir TMP

gatk --java-options "-Xmx${task.memory.giga}g -Djava.io.tmpdir=TMP" GenotypeConcordance \
	--TRUTH_VCF "${row.vcf1}" \
	--TRUTH_SAMPLE "${row.sample}" \
	--CALL_VCF "${row.vcf2}" \
	--CALL_SAMPLE "${row.sample}_2" \
	--O "${meta.prefix?:""}${row.sample}"


cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">GT Concordance</entry>
        <entry key="sample">${row.sample}</entry>
        <entry key="vcf1">${row.vcf1}</entry>
        <entry key="vcf2">${row.vcf2}</entry>
        <entry key="versions">${getVersionCmd("gatk")}</entry>
</properties>
EOF
"""
}

process PLOT {
tag "${files.name}"
input:
	val(meta)
	path(files)
output:
	path("${meta.prefix?:""}concordance.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
"""
hostname 1>&2
${moduleLoad("r")}


xargs -L 50 -a "${files}" cat | grep '^SNP' | cut -f 13 | grep -v -F -w '?' | cat > snps.txt
xargs -L 50 -a "${files}" cat | grep '^INDEL' | cut -f 13 | grep -v -F -w '?' | cat > indels.txt

cat << EOF | R --vanilla

T1<- read.table("snps.txt", header = FALSE, sep = "\t",col.names=c('C'), colClasses = c("numeric"),stringsAsFactors=FALSE)
T2<- read.table("indels.txt", header = FALSE, sep = "\t",col.names=c('C'), colClasses = c("numeric"),stringsAsFactors=FALSE)

pdf("${meta.prefix?:""}concordance.pdf")
boxplot(T1\\\$C,T2\\\$C,main="${params.prefix} Condordance",sub="",xlab="Variant Type",ylab="Concordance",col=c("blue","red"),names=c("SNPs","INDELs"))
dev.off()
EOF


cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">GT Concordance</entry>
        <entry key="versions">${getVersionCmd("R")}</entry>
</properties>
EOF
"""
}

process ZIP {
executor "local"
input:
	val(meta)
	path(concordances)
	path(version)
	path(html)
	path(pdf)
output:
	path("${meta.prefix?:""}archive.zip"),emit:zip
script:
"""
cp "${concordances}" jeter.list
echo "${version}" >> jeter.list
echo "${html}" >> jeter.list
echo "${html}" >> jeter.list
echo "${pdf}" >> jeter.list

zip -9 -@ -j "${meta.prefix?:""}archive.zip" < jeter.list
rm jeter.list
"""
}
