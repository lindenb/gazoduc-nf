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


include {runOnComplete;moduleLoad;dumpParams} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {MULTIQC} from '../../subworkflows/multiqc/multiqc.nf'


if( params.help ) {
    dumpParams(params);
    exit 0
}  else {
    dumpParams(params);
}




workflow {
	ch = AF_FREQ_GNOMAD_VS_VCF([:],
		params.genomeId,
		params.vcf,
		Channel.fromPath(params.bed),
		params.sample2pop
		)

	}


workflow AF_FREQ_GNOMAD_VS_VCF {
	take:
		meta
		genomeId
		vcf
		bed
		sample2pop
	main:

		version_ch = Channel.empty()

		interval_ch = bed.splitCsv(header: false,sep:'\t',strip:true).
			map{T->[(T.size() < 4 ? "ALL":T[3]), T[0]+"\t"+T[1]+"\t"+T[2] ]}.
			groupTuple()

		plot_ch = PLOT_INTERVAL([:], genomeId, vcf, sample2pop, interval_ch)
		version_ch = version_ch.mix(plot_ch.version)

		MULTIQC(plot_ch.output.collect())

		version_ch = MERGE_VERSION("vsgnomad",version_ch.collect())
	emit:
		version = version_ch
	}

runOnComplete(workflow)


process PLOT_INTERVAL {
tag "${title} n=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	val(vcf)
	val(sample2pop)
	tuple val(title),val(L)
output:
	path("${params.prefix?:""}${title}.*"),emit:output
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
	def reference = genome.fasta
	def max_diff=0.1
	def gnomad = genome.gnomad_genome
	def description = "compare VCF allele frequencies with <code>${params.af_tag}</code> in gnomad."
"""
hostname 1>&2
${moduleLoad("bcftools bedtools jvarkit")}
set -o pipefail
mkdir -p TMP
set -x

# merge exon bed
cat << EOF |  cut -f1,2,3 | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/intervals.bed
${L.join("\n")}
EOF

# get gnomad header
bcftools view --header-only "${gnomad}" > TMP/header.gnomad.vcf

# check INFO/af_xx is present in header
grep -w -F '${params.af_tag}' TMP/header.gnomad.vcf 1>&2

# convert bed for gnomad
java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.gnomad.vcf --column 1 --convert SKIP TMP/intervals.bed > TMP/gnomad.bed

# convert the bed for the vcf
bcftools view --header-only "${vcf}" > TMP/header.user.vcf
java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.user.vcf --column 1 --convert SKIP TMP/intervals.bed > TMP/user.bed


# extract BAD intervals for gnomad
bcftools view -O u -G ${params.exclude_gnomad_expression} \
		--regions-file TMP/gnomad.bed \
		"${gnomad}" |\
	bcftools query -f '%CHROM\t%POS0\t%END\\n' > TMP/exclude.01.bed

# also multi allelic are split in gnomad
bcftools query -f '%CHROM\t%POS0\t%REF\\n' --regions-file TMP/gnomad.bed "${gnomad}" |\
	sort -T TMP | uniq -d |\
	awk '{printf("%s\t%s\t%d\\n",\$1,\$2,int(\$2)+1);}' >> TMP/exclude.01.bed


# extract BAD intervals for vcf
bcftools view -O u -G ${params.exclude_vcf_expression} \
		--regions-file TMP/user.bed \
		"${vcf}" |\
	bcftools query -f '%CHROM\t%POS0\t%END\\n' > TMP/exclude.02.bed

# merge BAD intervals for gnomad
cat TMP/exclude.01.bed TMP/exclude.02.bed |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.gnomad.vcf  --column 1 --convert SKIP > TMP/exclude.gnomad.bed

# prevent empty file
if test ! -s TMP/exclude.gnomad.bed ; then
	awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' '${reference}.fai' | tail -n 1  > TMP/exclude.gnomad.bed
fi


# merge BAD intervals for vcf
cat TMP/exclude.01.bed TMP/exclude.02.bed |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.user.vcf  --column 1 --convert SKIP > TMP/exclude.user.bed

# prevent empty file
if test ! -s TMP/exclude.user.bed ; then
	awk -F '\t' '{printf("%s\t0\t1\\n",\$1);}' '${reference}.fai' | tail -n 1  > TMP/exclude.user.bed
fi



# extract gnomad
bcftools view -O u -m 2 -M 2 --types '${params.variant_types}' \
		${params.gnomadViewOpt?:""} \
		--regions-file TMP/gnomad.bed \
		"${gnomad}" |\
	bcftools view -O u --targets-file ^TMP/exclude.gnomad.bed --targets-overlap 1 |\
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT,%INFO/${params.af_tag}\\n' |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.user.vcf  --column 1 --convert SKIP |\
	sort -T . -t ',' -k1,1 > TMP/gnomad.af.csv

# extract samples in vcf
bcftools query --list-samples "${vcf}" | sort | uniq > TMP/samples.vcf.txt

touch TMP/labels.tsv


# get a file sample/population
if ${file(sample2pop).name.equals("NO_FILE")} ; then
	bcftools query -l '${vcf}' | awk '{printf("%s\tALL\\n",\$1);}' > TMP/samples.txt
else
	cp -v "${sample2pop}" TMP/samples.txt
fi


cut -f2 TMP/samples.txt  | uniq | sort | uniq | while read POP
do
	# subset samples for this population
	awk -v P=\${POP} -F '\t' '(\$2==P) {print \$1;}' TMP/samples.txt |\
		sort | uniq | comm -12 TMP/samples.vcf.txt - > TMP/pop.samples.txt 
	
	if ! [ -s TMP/pop.samples.txt ] ; then
		continue
	fi

	# extract samples
	bcftools norm -f "${reference}" --regions-file TMP/user.bed --multiallelics -any -O u "${vcf}" |\\
	bcftools view -O u -m 2 -M 2 \\
		--types '${params.variant_types}' \\
		${params.userViewOpt?:""} \\
		--samples-file TMP/pop.samples.txt | \\
	bcftools view -i 'ALT!="*"' -O u --targets-file ^TMP/exclude.user.bed --targets-overlap 1 |\
	bcftools +fill-tags -O u -- -t AF |\
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT,%INFO/AF\\n' |\
	sort -T TMP -t ',' -k1,1 > TMP/user.af.csv
	
	join -o '1.1,2.1,1.2,2.2' -t ',' -1 1 -2 1 TMP/gnomad.af.csv TMP/user.af.csv > TMP/join.csv

	# unmatched lines output is VAR1,VAR2,AF1,AF2
	join -v 1 -t ',' -1 1 -2 1 TMP/gnomad.af.csv TMP/user.af.csv | awk -F, '{printf("%s,NA,%s,0.0\\n",\$1,\$2);}'  >> TMP/join.csv
	join -v 2 -t ',' -1 1 -2 1 TMP/gnomad.af.csv TMP/user.af.csv | awk -F, '{printf("NA,%s,0.0,%s\\n",\$1,\$2);}'  >> TMP/join.csv
	tail TMP/join.csv 1>&2

	# remove '.' from data, sometimes shit happens
	awk -F, '!(\$3=="." || \$4==".")' TMP/join.csv > TMP/join.csv.bis
	mv -v TMP/join.csv.bis TMP/join.csv


	# high differences
	awk -vPOP=\${POP} -F, '{T=\$1;if(T=="NA") T=\$2; split(T,a,/[\t]/); W=sprintf("%s:%d:%s",a[1],a[2],POP);A=\$3*1.0;B=\$4*1.0;D=A-B; if(D<0) D=D*-1;if(D>=${max_diff}) {printf("%f\t%f\t%s\\n",\$3,\$4,W);}  }' TMP/join.csv |\
		sort | uniq >> TMP/labels.tsv
	

	# publish if not empty
	if [ -s "TMP/join.csv" ] ; then

		cut -d, -f3,4  "TMP/join.csv" > "TMP/\${POP}.join.csv"

		echo -n "\${POP}\tTMP/\${POP}.join.csv" >> TMP/manifest.txt

		# title for pdf with number of samples
		awk -v P=\${POP} 'END{printf("\t%s (N=%d)\\n",P,NR);}' TMP/pop.samples.txt >> TMP/manifest.txt
	fi

	rm TMP/pop.samples.txt TMP/join.csv
done

sort -T TMP TMP/labels.tsv | uniq > TMP/labels.tsv.bis && mv TMP/labels.tsv.bis TMP/labels.tsv

cat << __EOF__ > jeter.R
TT<-read.table("TMP/manifest.txt",header = FALSE,sep="\t",comment.char="",col.names=c("title","path","title2"),stringsAsFactors=FALSE)
head(TT)

labels <-read.table("TMP/labels.tsv",header = FALSE,sep="\t",comment.char="",col.names=c("x","y","label"),stringsAsFactors=FALSE,colClasses=c("numeric","numeric","character"))
head(labels)

colors <-  rainbow(nrow(TT),alpha=0.6)

maxXY <- 0.0

# get max AF in all files
for(i in c(1:nrow(TT))) {
        T2 <- read.table(TT[i,2],header = FALSE,sep=",",comment.char="",col.names=c("X1","X2"),colClasses=c("numeric","numeric"))
        if(nrow(T2)==0) continue;
	maxXY = max(T2[,1])
	maxXY = max(T2[,2])
	}

# extend a little so see max points
maxXY <- maxXY + maxXY*0.01

if(maxXY==0.0 || maxXY > 0.8) {
	maxXY <- 1.001
	}

# max af was specified by user
if(${params.max_af} > 0) {
	maxXY = ${params.max_af}
	}

png("TMP/out.png")

for(i in c(1:nrow(TT))) {
	T2<-read.table(TT[i,2],header = FALSE,sep=",",comment.char="",col.names=c("X1","X2"),colClasses=c("numeric","numeric"))
	if(nrow(T2)==0) continue;
	T2<-as.matrix(T2)
	if(i==1) {
		plot(T2,
			type = "p",
			main="${title}",
			sub ="INFO/${params.af_tag}",
			xlab="AF ${file(gnomad).name}",
			ylab="AF ${file(vcf).name}",
			las=2,
			xlim=c(0,maxXY),
			ylim=c(0,maxXY),
			col= colors[i],
			pch=16
			)
		# diagonal
		segments(0.,0,1.0,1.0,col="gray",lty="dotted")

		# plot vertical and horizontal delim for rares frequencies
		rareAF <- 0.01
		abline(v = rareAF, col="gray",lty="dotted")
		abline(h = rareAF, col="gray",lty="dotted")

		# labels
		if(nrow(labels)>0) {
			text(x=labels[,1],y=labels[,2],labels=labels[,3], cex=0.4)
			}
	} else {
		par(new = TRUE)
		plot(T2,
			type = "p",
			axes=FALSE,
			ann=FALSE,
			col= colors[i],
			xlim=c(0,maxXY),
			ylim=c(0,maxXY),
			pch=16
			)
		
	}

}

legend("bottomright",legend=TT[,3],title="population",pch=16,col=colors) 
dev.off()
__EOF__


R --vanilla --quiet < jeter.R



mv -v "TMP/out.png" "${params.prefix?:""}${title}.png"

cut -f 3 TMP/labels.tsv | sort | uniq > "${params.prefix?:""}${title}.differences.intervals"

# extract highest AF, sort
awk -F '\t' '{M=(\$1*1.0 - \$2*1.0);if(M<0) {M=-M;}  printf("%f\t%s\\n",M,\$0);}' TMP/labels.tsv | sort -T TMP -t '\t' -k1,1gr | cut -f2- | head -n 10  > TMP/jeter1.txt


if test -s TMP/jeter1.txt ; then


cat << __EOF__ > TMP/jeter.html
<!--
parent_id: af_section_diff
parent_name: "VCF vs Gnomad. Major AF differences"
parent_description: "A few points for each gene where the AF in the VCF was different from the AF in Gnomad"
id: '${title}_table_diff'
section_name: '${title} differences'
description: 'Some of the worst differences (&Delta;AF greater than <code>${max_diff}</code>) between the VCF and gnomad in <b>${title}</b>.'
-->
<pre>
__EOF__

awk '(NR==1) {print "GNOMAD\tUSER\tVARIANT";} {print}'  TMP/jeter1.txt | column -t >> TMP/jeter.html

echo "</pre>" >> TMP/jeter.html

mv TMP/jeter.html "${params.prefix?:""}${title}.table_mqc.html"

fi


##
## create MULTIQC CONFIG
##
cat << EOF > "${params.prefix?:""}${title}.multiqc_config.yaml"
custom_data:
  af_${title}:
    parent_id: af_section
    parent_name: "VCF vs Gnomad"
    parent_description: "${description}"
    section_name: "${title}"
    description: "Compare allele frequencies for <b>${title}</b> in <code>${vcf}</code> ${params.max_af>0?"max AF:<code>${params.max_af}</code>":""}."
sp:
  af_${title}:
    fn: "${params.prefix?:""}${title}.png"
ignore_images: false
EOF




##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">plot AF vs AF</entry>
</properties>
EOF
"""
}
