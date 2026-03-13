
process DIGEST_SAMPLESHEET {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1 ),path(vcf)
	tuple val(meta ),path(samplesheet)
output:
	tuple val(meta ),path("*.cases.txt"),emit:cases
	tuple val(meta ),path("*.controls.txt"),emit:controls
	tuple val(meta ),path("*.all.samples.txt"),emit:samples
	tuple val(meta ),path("*.plink.ped"),emit:plink_ped
	tuple val(meta ),path("*.sample2population.tsv"),emit:sample2population
	path("versions.yml"),emit:versions
script:
    def prefix = task.ext.id?:"${meta.id}"
"""
hostname 1>&2
mkdir -p TMP
set -x


# samples in VCF
bcftools query -l "${vcf}" | sort | uniq > TMP/in.vcf.samples.txt
test -s TMP/in.vcf.samples.txt


cat << 'EOF' > TMP/jeter.R

SN= read.table("TMP/in.vcf.samples.txt",header=TRUE,sep="\\t", stringsAsFactors=FALSE,col.names=c("sample"))


T1 <- read.table("${samplesheet}",header=TRUE,sep="\\t", stringsAsFactors=FALSE)

required_columns <- c("sample", "sex", "status")
missing_columns <- setdiff(required_columns, colnames(T1))

if (length(missing_columns) > 0) {
  stop(paste("Erreur : Les colonnes suivantes sont manquantes dans le fichier :", paste(missing_columns, collapse=", ")))
}


if(!"population" %in% colnames(T1)) {
  # Copy the content of 'status' column into the new 'population' column
  T1\$population <- T1\$status
}
meta1
T1\$population[T1\$population == ""] <- "NA"
T1\$population[T1\$population == "."] <- "NA"



T1 <- T1[T1\$sample %in% SN\$sample,]

if (any(duplicated(T1\$sample))) {
  stop("DUPLICATE SAMPLES")
}

T2 <- T1[T1\$status=='case',]\$sample
write.table(T2,"TMP/cases.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)
T2 <- T1[T1\$status=='control',]\$sample
write.table(T2,"TMP/controls.txt",sep="\t", row.names=FALSE, quote=FALSE, col.names = FALSE)


# Conversion des valeurs de sex pour PLINK (1 = male, 2 = female, 0 = unknown)
T1\$sex <- ifelse(T1\$sex == "male", 1, ifelse(T1\$sex == "female", 2, 0))
T1\$status <- ifelse(T1\$status == "case" , "1", ifelse(T1\$status == "control" , "0" , "NA"))

# header
T1\$FID <- T1\$sample
T1\$IID <- T1\$sample


T2<-T1[,c("sample","population")]
write.table(T2,file="TMP/sample2pop.txt", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

T1 <- T1[, c("FID", "IID", "sex", "${params.status}")]
write.table(T1,file="TMP/plink.ped", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
EOF

R --no-save --vanilla < TMP/jeter.R


sort TMP/cases.txt > ${prefix}.cases.txt
sort TMP/controls.txt > ${prefix}.controls.txt
cat TMP/cases.txt TMP/controls.txt > ${prefix}.all.samples.txt
mv TMP/plink.ped ${prefix}.plink.ped
mv TMP/sample2pop.txt ${prefix}.sample2population.tsv
touch versions.yml
"""
stub:
    def prefix = task.ext.id?:"${meta.id}"
"""
touch versions.yml  ${prefix}.cases.txt ${prefix}.controls.txt ${prefix}.all.samples.txt ${prefix}.plink.ped ${prefix}.sample2population.tsv
"""
}

/****************************************************************************************************/


process MDS_TO_COVARIATES {
executor "local"
input:
	tuple val(meta),path(mds)
output:
	tuple val(meta),path("*.tsv"),emit:covariates
	path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id}.mds2covar"
"""
awk 'BEGIN{printf("FID\tIID\tY1\tY2\tY3\\n");} (NR>1) {printf("%s\t%s\t%s\t%s\t%s\\n",\$2,\$2,\$4,\$5,\$6);}' "${mds}" > "${prefix}.covariates.tsv"
touch versions.yml
"""
stub:
    def prefix = task.ext.prefix?:"${meta.id}.mds2covar"
"""
touch versions.yml ${prefix}.tsv
"""
}

/****************************************************************************************************/

process UPDATE_PGEN {
label "process_single"
tag "${meta.id}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(plink_pedigree)
	tuple val(meta ),path(pgen),path(psam),path(pvar)
output:
    tuple val(meta),path("*.pgen"),path("*.psam"),path("*.pvar"),emit:pgen
	path("versions.yml"),emit:versions
script:
	def prefix = task.prefix?:"${meta.id}.update"
	def status = task.ext.status?:"status" 
"""
mkdir -p TMP
#
# UPDATE file psam with the information of the pedigree
#

# sort on name, keep original order in 1st column: INDEX/SAMPLE
awk -F '\t' '(NR>1) {printf("%d\t%s\\n",NR,\$1);}' '${psam}' |\\
	sort -T TMP -t '\t' -k2,2 > TMP/jeter.a

head TMP/jeter.a 1>&2

#
# sort pedigree SAMPLE,FATHER,MOTHER
#
cut -f2,3,4 "${plink_pedigree}" |\\
	tail -n +2 |\\
	sort -T TMP -t '\t' -k1,1 > TMP/jeter.b

head TMP/jeter.b 1>&2

#
# JOIN, sort on ORIGINAL ORDER
#
join -t '\t' -1 2 -2 1 -e 'NA' -a 1 -o '1.1,1.2,2.2,2.3' TMP/jeter.a TMP/jeter.b |\\
	sort -t '\t' -k1,1n |\\
	cut -f 2- |\\
	awk -F '\t' 'BEGIN{printf("#FID\tIID\tSEX\t${status}\\n");} {printf("%s\t%s\t%s\t%s\\n",\$1,\$1,\$2,\$3);}' > TMP/new.psam

ln -s ${pgen} ${prefix}.pgen
ln -s ${pvar} ${prefix}.pvar
mv TMP/new.psam ${prefix}.psam


touch versions.yml
"""

stub:
	def prefix = task.prefix?:"${meta.id}.update"
"""
touch versions.yml ${prefix}.pgen ${prefix}.pvar ${prefix}.psam
"""
}


/****************************************************************************************************/

process FUNCTIONAL_ANNOTATION_SCORES {
executor "local"
input:
	tuple val(meta),path(user_custom_scores)
output:
	tuple val(meta),path("scores.tsv"),emit:tsv
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"scores"
"""
test -s ${user_custom_scores}
tr -s " " < "${user_custom_scores}" | tr " " "\t" > jeter.tsv

awk -F '\t' '(NF!=3 || \$1=="" || \$2=="" || \$3=="")'  jeter.tsv > valid.txt
# error if any column 1-3 is empty
test ! -s valid.txt
rm valid.txt

mv jeter.tsv ${prefix}.tsv
touch versions.yml
"""

stub:
def prefix = task.ext.prefix?:"scores"
"""
touch ${prefix}.tsv versions.yml
"""
}


/****************************************************************************************************/

process MAKE_FUNCTIONAL_ANNOT_PER_CTG {
tag "chr${contig}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"   
afterScript "rm -rf TMP"
input:
	tuple val(meta1),path(annotations)
	tuple val(meta2),path(gencode),path(gencode_tbi)
    tuple val(contig),path(vcf),path(vcf_tbi)
output:
    tuple val("functional"),val(contig),path("OUT/manifest.tsv"),emit:output
	path("versions.yml"),emit:versions
when:
	!(params.skip_XY && contig.matches("(chr)?[XY]"))
script:
	def jvm =tak.ext.jvm?:"-Djava.io.tmpdir=TMP "
	def jvarkit = task.ext.jvarkit?:"java -jar  ${jvm} \${HOME}/jvarkit.jar"
	def contig = meta.contig
	def freq = task.ext.freq?:""
	verify(!isBlank(contig),"${task.process} : contig is blank")
"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT

bcftools view -O v '${vcf}' |\\
	${jvarkit} regeniefunctionalannot \\
		--annotations "${annotations}" \\
		--kg '${gencode}' \\
		-f "${freq}" |\\
	${jvarkit} regeniemakeannot \\
		-m "${annotations}" \\
		--prefix "chr${contig}chunk" \\
		--reserve 20 \\
		-o \${PWD}/OUT \\
		--gzip \\
		-N 5000

touch versions.yml
"""
stub:
"""
touch versions.yml
"""
}

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/

process MAKE_BED {
tag "chr${contig} ${select_bed.name}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
input:
	tuple val(meta),path(select_bed)
    tuple val(contig),path(vcf),path(tbi)
output:
        tuple val("user_bed"),val(contig),path("OUT/manifest.tsv"),emit:output
		path("versions.yml"),emit:versions
when:
	def contig = "${meta.contig}"
	!(params.skip_XY && contig.matches("(chr)?[XY]"))
script:
	def jvm =tak.ext.jvm?:"-Djava.io.tmpdir=TMP "
	def jvarkit = task.ext.jvarkit?:"java -jar  ${jvm} \${HOME}/jvarkit.jar"
	def min_length= task.ext.min_bed_length?:"0" //(params.min_bed_length?:0)
"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT

# extract DICT only
bcftools view --header-only  -O z -o TMP/dict.vcf.gz '${vcf}'

# rename contig in user bed
${select_bed.name.endsWith(".gz")?"gunzip -c ":"cat"} "${select_bed}" |\\
	${jvarkit} bedrenamechr -R TMP/dict.vcf.gz --column 1 > TMP/jeter.bed

# prevent empty file
if test ! -s TMP/jeter.bed
then
	echo -e "XXXXX\t0\t1" > TMP/jeter.bed
fi

bcftools view --regions-file TMP/jeter.bed -O v '${vcf}' |\\
	${jvarkit} regeniebedannot \\
		--bed "${select_bed}" \\
		--min-length ${min_length} \\
		${params.skip_XY?"--noXY":""} \\
		-f ${params.freq} |\\
	${jvarkit} regeniemakeannot \\
		--prefix "chr${contig}_bed_chunk" \\
		-o \${PWD}/OUT \\
		--reserve 10 \\
		--gzip \\
		-N 5000

touch versions.yml
"""
stub:

"""
touch versions.yml OUT/manifest.tsv
"""
}

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/

process MAKE_SLIDING {
tag "chr${contig} ${win_size}/${win_shift}"
label "process_single"
conda "${moduleDir}/../../conda/bioinfo.01.yml"   
afterScript "rm -rf TMP"
input:
        tuple val(contig),path(vcf_files),val(win_size),val(win_shift)
output:
        tuple val("sliding_${win_size}_${win_shift}"),val(contig),path("OUT/manifest.tsv"),emit:output
when:
	!(params.skip_XY && contig.matches("(chr)?[XY]"))
script:
        def vcf = vcf_files.find{it.name.endsWith(".bcf") || it.name.endsWith(".vcf.gz")}

"""
set -o pipefail
mkdir -p TMP
mkdir -p OUT

bcftools view -O v '${vcf}' |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regenieslidingannot \\
		--window-size "${win_size}" \\
		--window-shift "${win_shift}" \\
		-f ${params.freq} |\\
	java -Djava.io.tmpdir=TMP -jar "\${HOME}/packages/jvarkit/dist/jvarkit.jar" regeniemakeannot \\
		--prefix "chr${contig}_${win_size}_${win_shift}_chunk" \\
		-o \${PWD}/OUT \\
		--reserve 20 \\
		--gzip \\
		-N 5000
	
"""
}




process MERGE_AND_PLOT {
tag "${title}"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
afterScript "rm -rf TMP"
label "process_single"
input:
	tuple val(meta1),path(fai)
	tuple val(meta1),path(dict)
	tuple val(meta ),val(L) //path("INPUT/*")
output:
	tuple val(meta),path("*.results.tsv.gz"),emit:tsv
	tuple val(meta),path("*.png"),emit:images
	tuple val(meta),path("*.report.txt"),emit:ascii
script:
	def title = task.ext.title?:"${meta.id}"
	def prefix = task.ext.prefix?:"${meta.id}"
	def regex_contig= "0-9X"
	def hline1 = 5
	def hline2 = hline1 + 1
"""
mkdir -p TMP
set -x

JD1=`which jvarkit`
echo "JD1=\${JD1}" 1>&2
# directory of jvarkit
JD2=`dirname "\${JD1}"`
# find the jar itself
JVARKIT_JAR=`find "\${JD2}/../.." -type f -name "jvarkit.jar" -print -quit `

cp -v "${moduleDir}/Minikit2.java" TMP/Minikit.java
javac -sourcepath TMP -cp "\${JVARKIT_JAR}" -d TMP TMP/Minikit.java

# do not use chrY because regenie merge it with chrX (see doc)
cut -f1,2 '${fai}' | awk -F '\t' '( \$1 ~ /^(chr)?[${regex_contig}]+\$/ )' > TMP/jeter.fai


cat << EOF > TMP/jeter.list
${L.join("\n")}
EOF

##cat INPUT/*.regenie.gz 

xargs -a TMP/jeter.list -L 1 gunzip -c |\\
	grep  '^CHROM' -m1 > TMP/jeter.txt

test -s TMP/jeter.txt

set -o pipefail

#cat INPUT/*.regenie.gz |\\
#	gunzip -c


xargs -a TMP/jeter.list -L 1 gunzip -c |\\
	grep -v "^#" |\\
	grep -v '^CHROM' |\\
	LC_ALL=C sort -T TMP  --buffer-size=${task.memory.mega}M  -t ' ' -k1,1V -k2,2n |\\
	uniq >> TMP/jeter.txt

cat TMP/jeter.txt | tr " " "\t" | gzip --best > TMP/${title}.results.tsv.gz


java -cp "\${JVARKIT_JAR}:TMP" Minikit -R "${fasta}" -o TMP < TMP/jeter.txt
rm TMP/jeter.txt

find ./TMP -type f 1>&2

cat << '__EOF__' > TMP/jeter.R
mft <- read.table("TMP/manifest.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE, colClasses=c("character","character","character"),)
head(mft)

fai <- read.table("TMP/jeter.fai", header=FALSE, sep="\t", stringsAsFactors=FALSE, colClasses=c("character","double"))
colnames(fai) <- c("contig", "size")

pos2index <- function(CHROM, POS) {  
  total_size <- 0.0
  for (i in 1:nrow(fai)) {
    if (fai\$contig[i] == CHROM) {
      return(total_size + POS)
    }
    total_size <- total_size + fai\$size[i]
  }
    stop("Chromosome non trouve")
}



chrom_positions <- c(0, cumsum(fai\$size))
chrom_colors <- rep(c("azure", "azure2"), length.out=nrow(fai))


last_chrom <- fai\$contig[nrow(fai)]
last_size <- fai\$size[nrow(fai)]
max_genome_size <- pos2index(last_chrom, last_size)


for(i in seq_len(nrow(mft))) {

data <- read.table(mft[i,]\$filename, header=TRUE, sep=" ", stringsAsFactors=FALSE)
data\$LOG10P <-  as.double(data\$LOG10P)

# save for html report
data2 <- data[is.finite(data\$LOG10P) & as.double(data\$LOG10P) >= ${hline2} , ]
fileout <- paste("TMP/${title}.freq_",mft[i,]\$freq,"_",mft[i,]\$test_name,"_",mft[i,]\$annot,".regenie.report.txt",sep="")
write.table(data2,fileout,sep="\t", row.names=FALSE, quote=FALSE, col.names = TRUE)



data\$x <- mapply(pos2index, data\$CHROM, data\$GENPOS)

fileout <- paste("TMP/${title}.freq_",mft[i,]\$freq,"_",mft[i,]\$test_name,"_",mft[i,]\$annot,".regenie.manhattan.png",sep="")

png(fileout, width=1500, height=500,units="px")

x_lim <-  c(0,max_genome_size)
y_lim <- c(0,1+max(data\$LOG10P))

plot(	NULL,
	main=paste("Test:",mft[i,]\$test_name," AF:",mft[i,]\$freq," Annot:",mft[i,]\$annot,sep=""),
	sub= "${params.prefix}regenie",
	xlab="${fasta}",
	ylab="-log10(PVALUE)",
	xlim = x_lim,
	ylim = y_lim
	)

for (k in seq_along(fai\$contig)) {
  rect(chrom_positions[k], par("usr")[3], chrom_positions[k+1], par("usr")[4], col=chrom_colors[k], border=NA)
}

# real plot
points(data\$x,
	data\$LOG10P,
	pch=20,
	col="black",
	xlim = x_lim,
	ylim = y_lim
	)

# vertical bars
abline(v=chrom_positions, col="black", lty=2)
abline(h=${hline2},col="red",lty=2)
abline(h=${hline1},col="green",lty=2)

dev.off()

pvector <- data\$LOG10P
# I don't understand why running pow(10,x) and then log10(y) gives a good chart. I hate R.
pvector <- 10^(-pvector)
pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector>0]
if(length(pvector) > 0) {
fileout <- paste("TMP/${title}.freq_",mft[i,]\$freq,"_",mft[i,]\$test_name,"_",mft[i,]\$annot,".regenie.qqplot.png",sep="")
png(fileout, width=500, height=500,units="px")


# Observed and expected
o = -log10(sort(pvector,decreasing=FALSE))
e = -log10( ppoints(length(pvector) ))

plot(x=e,
	y=o,
	pch=20, 
	xlim=c(0, max(e)),
	ylim=c(0, max(o)), 
        xlab=expression(Expected~~-log[10](italic(p))), 
        ylab=expression(Observed~~-log[10](italic(p))),
	main=paste("Test:",mft[i,]\$test_name," AF:",mft[i,]\$freq," Annot:",mft[i,]\$annot,sep=""),
	sub= "${params.prefix}regenie"
	)

abline(0,1,col="red")
dev.off()
}
}
__EOF__

R --no-save < TMP/jeter.R

# normalize columns
find TMP -type f -name "*.regenie.report.txt" | while read F
do
	column -t --separator \$'\\t' "\${F}" > TMP/jeter.column
	mv TMP/jeter.column "\${F}"
done

find TMP/ -type f 1>&2
mv TMP/${prefix}.results.tsv.gz ./
mv TMP/*.png ./
mv TMP/*.regenie.report.txt ./
touch versions.yml
"""
stub:
"""
touch versions.yml
"""
}