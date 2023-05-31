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

def gazoduc = gazoduc.Gazoduc.getInstance(params).putDefaults().putReference()


gazoduc.build("vcf","NO_FILE").
        desc("Input VCF").
        required().
        existingFile().
        put()

gazoduc.build("sample2pop","NO_FILE").
        desc("optional tab delimited file: SAMPLE(tab)POPULATION").
        put()

gazoduc.build("af_tag","AF_popmax").
        desc("INFO/TAG in gnomad").
        put()

gazoduc.build("bed","NO_FILE").
        desc("bed file CHROM/START/END/TITLE . Results will be grouped by title").
	required().
	existingFile().
        put()


params.gnomadViewOpt=""
params.userViewOpt=""


include {runOnComplete;moduleLoad} from '../../modules/utils/functions.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.nf'
include {SIMPLE_PUBLISH_01} from '../../modules/utils/publish.simple.01.nf'

if( params.help ) {
    gazoduc.usage().
	name("plot AF vcf vs AF gnomad").
	desc("Compare AF in vcf vs gnomad").
	print();
    exit 0
} else {
   gazoduc.validate();
}



workflow {
	ch = AF_FREQ_GNOMAD_VS_VCF( params, params.reference, params.gnomad, params.vcf, Channel.fromPath(params.bed), params.sample2pop )

	SIMPLE_PUBLISH_01(params, Channel.empty().mix(ch.zip).mix(ch.version).mix(ch.pdf).collect())
	}


workflow AF_FREQ_GNOMAD_VS_VCF {
	take:
		meta
		reference
		gnomad
		vcf
		bed
		sample2pop
	main:

		version_ch = Channel.empty()

		interval_ch = bed.splitCsv(header: false,sep:'\t',strip:true).
			map{T->[(T.size() < 4 ? "ALL":T[3]), T[0]+"\t"+T[1]+"\t"+T[2] ]}.
			groupTuple()

		plot_ch = PLOT_INTERVAL(meta, reference, gnomad, vcf, sample2pop, interval_ch)
		version_ch = version_ch.mix(plot_ch.version)

		zip_ch = ZIP_IT(meta, plot_ch.output.collect(), plot_ch.diff.collect())
		version_ch = version_ch.mix(zip_ch.version)

		version_ch = MERGE_VERSION(meta, "vsgnomad", "vs gnomad",version_ch.collect())

	emit:
		pdf = zip_ch.pdf
		zip = zip_ch.zip
		version = version_ch
	}

runOnComplete(workflow)


process PLOT_INTERVAL {
tag "${title} n=${L.size()}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(reference)
	val(gnomad)
	val(vcf)
	val(sample2pop)
	tuple val(title),val(L)
output:
	path("${meta.prefix?:""}${title}.pdf"),emit:output
	path("${meta.prefix}${title}.differences.intervals"),emit:diff
	path("version.xml"),emit:version
script:
	def max_diff=0.1
"""
hostname 1>&2
${moduleLoad("bcftools bedtools jvarkit")}
set -o pipefail
mkdir -p TMP

# merge exon bed
cat << EOF |  cut -f1,2,3 | sort -T TMP -t '\t' -k1,1 -k2,2n | bedtools merge > TMP/intervals.bed
${L.join("\n")}
EOF

# convert the bed for gnomad
bcftools view --header-only "${gnomad}" > TMP/header.gnomad.vcf
grep -w -F '${meta.af_tag}' TMP/header.gnomad.vcf
java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.gnomad.vcf --column 1 --convert SKIP TMP/intervals.bed > TMP/gnomad.bed

# convert the bed for the vcf
bcftools view --header-only "${vcf}" > TMP/header.user.vcf
java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.user.vcf --column 1 --convert SKIP TMP/intervals.bed > TMP/user.bed


# extract BAD intervals for gnomad
bcftools view -O u -G -i '(FILTER!="PASS" && FILTER!=".") || TYPE!~"snp" || N_ALT>1' \
		--regions-file TMP/gnomad.bed \
		"${gnomad}" |\
	bcftools query -f '%CHROM\t%POS0\t%END\\n' > TMP/exclude.01.bed

# also multi allelic are split in gnomad
bcftools query -f '%CHROM\t%POS0\t%REF\\n' --regions-file TMP/gnomad.bed "${gnomad}" |\
	sort -T TMP | uniq -d |\
	awk '{printf("%s\t%s\t%d\\n",\$1,\$2,int(\$2)+1);}' >> TMP/exclude.01.bed


# extract BAD intervals for vcf
bcftools view -O u -G -i '(FILTER!="PASS" && FILTER!=".") || TYPE!~"snp" || N_ALT>1' \
		--regions-file TMP/user.bed \
		"${vcf}" |\
	bcftools query -f '%CHROM\t%POS0\t%END\\n' > TMP/exclude.02.bed

# merge BAD intervals for gnomad
cat TMP/exclude.01.bed TMP/exclude.02.bed |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.gnomad.vcf  --column 1 --convert SKIP > TMP/exclude.gnomad.bed

# merge BAD intervals for vcf
cat TMP/exclude.01.bed TMP/exclude.02.bed |\
	java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f TMP/header.user.vcf  --column 1 --convert SKIP > TMP/exclude.user.bed


# extract gnomad
bcftools view -O u -m 2 -M 2 --types snps \
		${meta.gnomadViewOpt?:""} \
		--regions-file TMP/gnomad.bed \
		"${gnomad}" |\
	bcftools view -O u --targets-file ^TMP/exclude.gnomad.bed --targets-overlap 1 |\
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT,%INFO/${meta.af_tag}\\n' |\
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
	bcftools view -O u -m 2 -M 2 \
		--types snps \
		${meta.userViewOpt?:""} \
		--samples-file TMP/pop.samples.txt \
		--regions-file TMP/user.bed "${vcf}" |\
	bcftools view -O u --targets-file ^TMP/exclude.user.bed --targets-overlap 1 |\
	bcftools +fill-tags -O u -- -t AF |\
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT,%INFO/AF\\n' |\
	sort -T . -t ',' -k1,1 > TMP/user.af.csv
	
	join -o '1.1,2.1,1.2,2.2' -t ',' -1 1 -2 1 TMP/gnomad.af.csv TMP/user.af.csv > TMP/join.csv

	# unmatched lines
	join -v 1 -t ',' -1 1 -2 1 TMP/gnomad.af.csv TMP/user.af.csv | awk -F, '{printf("%s,NA,%s,0.0\\n",\$1,\$2);}'  >> TMP/join.csv
	join -v 2 -t ',' -1 1 -2 1 TMP/gnomad.af.csv TMP/user.af.csv | awk -F, '{printf("NA,%s,0.0,%s\\n",\$1,\$2);}'  >> TMP/join.csv
	
	# high differences
	awk -F, '{T=\$1;if(T=="NA") T=\$2; split(T,a,/[\t]/); W=sprintf("%s:%d",a[1],a[2]);A=\$3*1.0;B=\$4*1.0;D=A-B; if(D<0) D=D*-1;if(D>=${max_diff}) {printf("%f\t%f\t%s\\n",\$3,\$4,W);}  }' TMP/join.csv |\
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

sort TMP/labels.tsv | uniq > TMP/labels.tsv.bis && mv TMP/labels.tsv.bis TMP/labels.tsv

cat << __EOF__ > jeter.R
TT<-read.table("TMP/manifest.txt",header = FALSE,sep="\t",comment.char="",col.names=c("title","path","title2"),stringsAsFactors=FALSE)
head(TT)

labels <-read.table("TMP/labels.tsv",header = FALSE,sep="\t",comment.char="",col.names=c("x","y","label"),stringsAsFactors=FALSE,colClasses=c("numeric","numeric","character"))
head(labels)

colors <-  rainbow(nrow(TT),alpha=0.6)

pdf("TMP/out.pdf")

for(i in c(1:nrow(TT))) {
	T2<-read.table(TT[i,2],header = FALSE,sep=",",comment.char="",col.names=c("X1","X2"))
	if(nrow(T2)==0) continue;
	T2<-as.matrix(T2)
	if(i==1) {
		plot(T2,
			type = "p",
			main="${title}",
			sub ="INFO/${meta.af_tag}",
			xlab="${file(meta.gnomad).name}",
			ylab="${file(meta.vcf).name}",
			las=2,
			xlim=c(0,1.0),
			ylim=c(0,1.0),
			col= colors[i],
			pch=16
			)
		segments(0.,0,1.0,1.0,col="gray",lty="dotted")
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
			xlim=c(0,1.0),
			ylim=c(0,1.0),
			pch=16
			)
		
	}

}

legend("bottomright",legend=TT[,3],title="population",pch=16,col=colors) 
dev.off()
__EOF__


R --vanilla --quiet < jeter.R



mv "TMP/out.pdf" "${meta.prefix}${title}.pdf"

cut -f 3 TMP/labels.tsv | sort | uniq > "${meta.prefix}${title}.differences.intervals"


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">plot AF vs AF</entry>
</properties>
EOF
"""
}


process ZIP_IT {
tag "N=${L.size()}"
afterScript "rm -rf TMP"
executor "local"
input:
	val(meta)
	val(L)
	val(L2)

output:
	path("${meta.prefix}all.zip"),emit:zip
	path("${meta.prefix}all.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
"""
mkdir -p TMP

cat ${L2.join(" ")} | sort -T . |uniq > TMP/${meta.prefix}differences.intervals

gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=${meta.prefix}all.pdf ${L.join(" ")}

zip -0 -j ${meta.prefix}all.zip  ${meta.prefix}differences.intervals ${meta.prefix}all.pdf


##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">join pdfs</entry>
</properties>
EOF
"""
}


