/*

Copyright (c) 2024 Pierre Lindenbaum

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

include { moduleLoad } from '../../modules/utils/functions.nf'
include {BED_CLUSTER_01 } from '../../modules/jvarkit/jvarkit.bedcluster.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {BCFTOOLS_CONCAT_01} from '../bcftools/bcftools.concat.01.nf'
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
//include {BED_WITH_VARIANTS_01} from '../../modules/bcftools/bed.with.variants.01.nf'

workflow SPLICEAI_01 {
	take:
		meta
		genomeId /* path to reference */
		vcf /* one indexed vcf or file with suffix .list with the path to the vcfs */
		pedigree /* jvarkit ped or NO_FILE */
		bed /* restrict to that bed or NO_FILE */
	main:
		if(!params.containsKey("spliceai")) throw new IllegalArgumentException("params.spliceai undefined");
		if(!params.spliceai.containsKey("splice_distance")) throw new IllegalArgumentException("params.spliceai.splice_distance undefined");

		version_ch = Channel.empty()
		exons_ch = DOWNLOAD_EXONS(genomeId, bed)
		version_ch = version_ch.mix(exons_ch.version)

		clusters_ch = BED_CLUSTER_01(["bed_cluster_method":"--jobs 50"], genomeId, exons_ch.bed)
		version_ch = version_ch.mix(clusters_ch.version)
		
		each_bed = clusters_ch.output.splitText().map{T->file(T.trim())}

		spliceai_ch = APPLY_SPLICEAI(genomeId, pedigree, vcf.combine(each_bed))
		version_ch = version_ch.mix(spliceai_ch.version)

		file_list_ch = COLLECT_TO_FILE_01([:],spliceai_ch.vcf.collect())
		version_ch = version_ch.mix(file_list_ch.version)

		concat_ch = BCFTOOLS_CONCAT_01([:],file_list_ch.output, bed)
		version_ch = version_ch.mix(concat_ch.version)

		plot1_ch = PLOT_IT([:], genomeId, concat_ch.vcf)
		version_ch = version_ch.mix(plot1_ch.version)

                version_ch = MERGE_VERSION("SpliceAI", version_ch.collect())
	emit:
		version = version_ch
		vcf = concat_ch.vcf
		pdf = plot1_ch.pdf
		
	}

process DOWNLOAD_EXONS {
	afterScript "rm -rf TMP"
	input:
		val(genomeId)
		path(bed)
	output:
		path("exons.bed.gz"),emit:bed
		path("version.xml"),emit:version
	script:
		if(!params.containsKey("spliceai")) throw new IllegalArgumentException("params.spliceai undefined");
		if(!params.spliceai.containsKey("splice_distance")) throw new IllegalArgumentException("params.spliceai.splice_distance undefined");

		def genome = params.genomes[genomeId]
		def reference = genome.fasta
		def distance = (params.spliceai.splice_distance as int)
		def url = genome.known_genes_url
	"""
	hostname 1>&2
	${moduleLoad("jvarkit bedtools")}
	set -o pipefail

	mkdir -p TMP

	if ${!bed.name.equals("NO_FILE")} ; then
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP "${bed}" |\
			sort -T TMP  -t '\t' -k1,1 -k2,2n |\
                	bedtools merge > TMP/jeter2.bed
	fi


	# url dans le help de spliceai
	wget -O - "${url}" |\
		gunzip -c |\
		awk -F '\t' '{N=int(\$9);X=${distance};split(\$10,S,/[,]/);split(\$11,E,/[,]/);for(i=1;i<=N;i++) {V=int(S[i]);if(i>1) printf("%s\t%d\t%d\\n",\$3,(V<X?0:V-X),V+X);V=int(E[i]);if(i<N) printf("%s\t%d\t%d\\n",\$3,(V<X?0:V-X),V+X);} }' |\
		java -jar \${JVARKIT_DIST}/bedrenamechr.jar -f "${reference}" --column 1 --convert SKIP |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		bedtools merge |\
		sort -T TMP -t '\t' -k1,1 -k2,2n |\
		${bed.name.equals("NO_FILE")?"":"bedtools intersect -u -a - -b TMP/jeter2.bed |"} \
		gzip --best > TMP/exons.bed.gz

	mv TMP/exons.bed.gz ./


	##################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">extract splice interval as bed file</entry>
		<entry key="url"><url>${url}</url></entry>
	        <entry key="splice_distance">${distance}</entry>
	</properties>
	EOF
	"""
	}

process APPLY_SPLICEAI {
	tag "${bed.name}"
        // conda "${moduleDir}/../../conda/spliceai.yml" BROKEN on our cluster
	conda "${params.conda.envs_path}/SPLICEAI"
	afterScript "rm -rf TMP"
	cpus 2
	input:
		val(genomeId)
		path(pedigree)
		tuple path(vcf),path(bed)
	output:
		path("contig.bcf"),emit:vcf
		path("contig.bcf.csi"),emit:index
		path("version.xml"),emit:version
	script:
		
		def genome= params.genomes[genomeId]
		def reference = genome.fasta
		def db =genome.spliceai_annotation_type
		def distance = params.spliceai.splice_distance?:50
	"""
	hostname 1>&2
	${moduleLoad("bcftools jvarkit")}
	set -o pipefail

	export OMP_NUM_THREADS=${task.cpus}

	mkdir -p TMP
	export TMPDIR=\${PWD}/TMP
		
	if ${vcf.name.endsWith(".list")} ; then	
		bcftools concat --remove-duplicates  --allow-overlaps --file-list "${vcf.toRealPath()}" --regions-file "${bed}" -O b -o TMP/jeter.bcf
	else
		bcftools view --regions-file "${bed}" -O b -o TMP/jeter.bcf "${vcf.toRealPath()}"
	fi

	if ${!pedigree.name.equals("NO_FILE")} ; then
 

		grep -v '^#' "${pedigree}" | cut -f 2 | sort -T TMP | uniq > TMP/a
		bcftools query -l TMP/jeter.bcf  | sort -T TMP | uniq > TMP/b
		comm -12 TMP/a TMP/b > TMP/jeter.samples
		test -s TMP/jeter.samples

		awk -F '\t' '(\$6=="affected" || \$6=="case")'  "${pedigree}" | cut -f 2 | sort -T TMP | uniq > TMP/cases1.txt
		awk -F '\t' '(\$6=="unaffected" || \$6=="control")'  "${pedigree}" | cut -f 2  | sort -T TMP | uniq > TMP/controls1.txt
		comm -12 TMP/b TMP/cases1.txt > TMP/cases.txt
		comm -12 TMP/b TMP/controls1.txt > TMP/controls.txt

		bcftools view -O b -o TMP/jeter2.bcf TMP/jeter.bcf
		mv TMP/jeter2.bcf TMP/jeter.bcf

	fi

	
	bcftools norm -f "${reference}" --multiallelics -both -O u TMP/jeter.bcf |\
		bcftools view  -e 'ALT="*" ' -O b -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf

	bcftools view TMP/jeter.bcf |\
	        spliceai -R "${reference}"  -A "${db}" -D ${distance} |\
        	bcftools view  -i 'INFO/SpliceAI!=""' -O b -o TMP/jeter2.bcf
	mv TMP/jeter2.bcf TMP/jeter.bcf


	if test -s TMP/cases.txt && test -s TMP/controls.txt ; then
		
		bcftools +contrast \
		     -0 "TMP/controls.txt" \
                     -1 "TMP/cases.txt" \
		     -a PASSOC,FASSOC,NASSOC,NOVELAL,NOVELGT \
		     -O b -o TMP/jeter2.bcf TMP/jeter.bcf

		mv TMP/jeter2.bcf TMP/jeter.bcf

		bcftools view TMP/jeter.bcf |\
		java -jar  \${JVARKIT_DIST}/jvarkit.jar vcftrio --pedigree "${pedigree}" |\
		bcftools view -O b -o TMP/jeter2.bcf
		mv TMP/jeter2.bcf TMP/jeter.bcf

	fi


	
	bcftools index TMP/jeter.bcf

	mv TMP/jeter.bcf ./contig.bcf 
	mv TMP/jeter.bcf.csi ./contig.bcf.csi

	##################
	cat <<- EOF > version.xml
	<properties id="${task.process}">
		<entry key="name">${task.process}</entry>
		<entry key="description">apply spliceai</entry>
		<entry key="db">${db}</entry>
		<entry key="pedigree">${pedigree}</entry>
		<entry key="bed">${bed}</entry>
	        <entry key="splice_distance">${distance}</entry>
		<entry key="bcftools.version">\$( bcftools --version-only)</entry>
		<entry key="spliceai.version">\$(spliceai --help | grep Version)</entry>
	</properties>
	EOF
	"""
	}

process PLOT_IT {
tag "${vcf.name}"
afterScript "rm -rf TMP"
input:
	val(meta)
	val(genomeId)
	path(vcf)
output:
	path("${params.prefix?:""}spliceai.pdf"),emit:pdf
	path("version.xml"),emit:version
script:
	def genome = params.genomes[genomeId]
"""
hostname 1>&2
${moduleLoad("bcftools jvarkit R")}
set -o pipefail
mkdir -p TMP
bcftools view "${vcf}" |\
	java -jar \${JVARKIT_DIST}/jvarkit.jar vcfmulti2oneinfo --info SpliceAI |\
	bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SpliceAI\\n' |\
	grep -v -F '0.00|0.00|0.00|0.00' |\
	tr "|" "\t" |\
	awk -F'\t' '(\$1 ~ /^(chr)?[0-9XY]+\$/)' > TMP/jeter.tsv


cat << '__EOF__' > TMP/jeter.R
a1 <- c("ACCEPTOR","DONOR")
a2 <- c("GAIN","LOSS")
a3 <- c("SCORE","POS")
header <- c("CHROM","LEN","REF","ALT","ALT2","GENE")
for(x3 in a3) {
	for(x1 in a1) {
		for(x2 in a2) {
			header <- c(header,paste(x3,x1,x2,sep="_"))
			}
		}
	}

DICT <- read.table("${genome.fasta}.fai",header = FALSE,sep="\t",comment.char="",col.names=c("CHROM","LENGTH","X1","X2","X3"),stringsAsFactors=FALSE)

DICT<-DICT[grep("^(chr)?[0-9XY]+\$",DICT\$CHROM),]


T1 <- read.table("TMP/jeter.tsv",header = FALSE,sep="\t",comment.char="",col.names=header,stringsAsFactors=FALSE)

pos2index <- function(dict,contig,pos) {
	i <- 0
	n <- 1
	prev <- 0
	while(n <= nrow(dict) )  {
		if(contig == DICT\$CHROM[n] ) {
			return (prev + pos );
			}
		prev <- prev + as.integer(DICT\$LEN[n])
		n <- n + 1
		}
	stop(contig)
	}

width <- 1000000
genome_size <- pos2index(DICT,DICT\$CHROM[nrow(DICT)],DICT\$LEN[nrow(DICT)])


indexes <- apply(T1,1,function(row) {pos2index(DICT,row[1],as.numeric(row[2]))})

head(indexes)
length(indexes)
head(T1[,7+0])
length(T1[,7+0])


pdf("TMP/jeter.pdf",width=40,height=25)
par(mfrow = c(4, 1))
y_max <- 1.1
for(col in c(0:3)) {
	plot(x=indexes,
		y=T1[,7+col],
		type="p",
		main=header[7+col],
		sub=paste("${params.prefix?:""}",header[7+col],sep=""),
		ylim=c(0, y_max),
		xlim=c(0,genome_size),
		xlab="Genome Position (${genome.name})",
		ylab=header[7+col],
		xaxt = "n"
		)

	gt8 <- T1[,7+col]>0.8
	
	if(nrow(T1[gt8,]) > 0 ) {
	text(
		x=indexes[gt8],
		y=T1[gt8,7+col],
		labels = T1[gt8,]\$GENE,
		cex= 1.0,
		col="blue"
		)
	}
	n<-1
	prev_x <- 0
	while(n <= nrow(DICT) )  {
		x<- pos2index(DICT,DICT\$CHROM[n],DICT\$LEN[n])
		abline( v=x ,col="blue")
		text((prev_x + x)/2.0 ,y_max,DICT\$CHROM[n],col="gray")
		prev_x <- x
		n <- n + 1 
		}
	}
dev.off()
__EOF__

R --vanilla < TMP/jeter.R
mv TMP/jeter.pdf  "${params.prefix?:""}spliceai.pdf"
	
##################
cat << EOF > version.xml
<properties id="${task.process}">
	<entry key="name">${task.process}</entry>
	<entry key="description">plot spliceai as manhattan</entry>
</properties>
EOF
"""
}
