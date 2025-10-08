/*

Copyright (c) 2025 Pierre Lindenbaum

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

include {SOMALIER_DOWNLOAD_SITES} from '../../../modules/somalier/download.sites'
include {EXTRACT_BAM            } from '../../../modules/somalier/extract'
include {RELATE                 } from '../../../modules/somalier/relate'


workflow SOMALIER_BAMS {
	take:
		meta
		fasta//MUST BE PROVIDED AS Channel.of() [meta,fasta]
		fai
		dict
		bams_ch // sample,bam,bai
		pedigree // pedigree for somalier
		user_sites //file or no file   [meta,vcf,vcf_idx]
	main:
		version_ch = Channel.empty()
		user_sites = user_sites.filter{it!=null && it.size()>1 && it[1]}


		fasta2 = bams_ch.count()
				.filter{it>1}
				.combine(fasta)
				.map{[it[1],it[2]]}

		/* sites PROVIDED */
		sites_vcf =  bams_ch.count()
			.filter{it>1}
			.combine(user_sites)
			.map{[it[1],it[2],it[3]]}


			fasta3 = user_sites
				.count()
				.filter{it==0}
                                .combine(fasta2)
                                .map{[it[1],it[2]]}


			SOMALIER_DOWNLOAD_SITES(fasta3,fai,dict)
			version_ch = version_ch.mix(SOMALIER_DOWNLOAD_SITES.out.versions)
			sites_vcf= sites_vcf.mix(SOMALIER_DOWNLOAD_SITES.out.vcf)

		sites_vcf = sites_vcf.first()


		ch1 = bams_ch.branch {
				tuple5: it.size()==5
				tuple3: it.size()==3
				other:true
				}


		ch1.other.map{throw new IllegalArgumentException("bad input for somalier ${it} size=${it.size()}");}


		ch2 = ch1.tuple3.combine(fasta2)
			.map{ [it[0],it[1],it[2],it[4], fai[1] ] }
			.mix( ch1.tuple5)


		EXTRACT_BAM(sites_vcf, ch2)
		version_ch = version_ch.mix(EXTRACT_BAM.out.versions.first())
	
		RELATE(
			fasta2,fai,dict,
			EXTRACT_BAM.out.output.map{it[1]}.collect().map{[[id:"somalier"],it]},
			pedigree
			)
		version_ch = version_ch.mix(RELATE.out.versions)

		PLOT_MATRIX(RELATE.out.output.map{[it[0],it[1].find{v->v.name.endsWith(".pairs.tsv")}]})
		version_ch = version_ch.mix(PLOT_MATRIX.out.versions)

	emit:
		output = RELATE.out.output
		versions = version_ch
		zip = RELATE.out.zip
		multiqc = RELATE.out.multiqc

}




process PLOT_MATRIX {
label "process_quick"
afterScript "rm -rf TMP"
conda  "${moduleDir}/../../../conda/bioinfo.01.yml"
tag "${meta.id}"
input:
	tuple val(meta),path(pairs)
output:
	tuple val(meta),path("*.pdf"),emit:pdf
	path("versions.yml"),emit:versions
script:
	def prefix = task.ext.prefix?:"${meta.id}"
	def title = task.ext.title?:"Relatedness Heatmap"
"""
mkdir -p TMP

cut -f1,2,3 "${pairs}" | sed 's/^#//' |\
	awk -F '\t' '(\$3!="-nan" && \$3!="-inf")'  > TMP/jeter.tsv 

cat << '__EOF__' | R --no-save 1>&2

# Read the data, treating '-nan' as NA
df <- read.table("TMP/jeter.tsv", header=TRUE, sep="\t", na.strings=c("-nan","-inf"))

head(df)

# Create a matrix with sample_b as rows and sample_a as columns
mat <- with(df, tapply(relatedness, list(sample_b, sample_a), function(x) {
  if(length(x)==0) return(NA)
  # If there are repeated pairs, take the mean; adjust as needed
  mean(x, na.rm=TRUE)
}))

summary(mat)

# Set up color palette (blue for low, white for mid, red for high, grey for NA)
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)

#min_rel <-  min(mat, na.rm=TRUE)
min_rel <- 0.0
#max(1.0,max(mat, na.rm=TRUE)
max_rel <-  1;0


print(min_rel)
print(max_rel)

# Open a PDF device
pdf("${prefix}.pdf", width=20, height=20)

# Plot the heatmap
image(
  1:ncol(mat), 1:nrow(mat), t(mat), 
  col=my_colors,
  axes=FALSE,
  xlab="Sample A",
  ylab="Sample B",
  zlim=c(min_rel, max_rel),
  main="${title}"
)
# Add NA as grey
#na_points <- which(is.na(t(mat)), arr.ind=TRUE)
#if(length(na_points) > 0) {
#  points(na_points[,1], na_points[,2], pch=15, col="grey80", cex=1.5)
#}
# Add axis labels
axis(1, at=1:ncol(mat), labels=colnames(mat), las=2, cex.axis=0.8)
axis(2, at=1:nrow(mat), labels=rownames(mat), las=2, cex.axis=0.8)
box()

dev.off()
__EOF__

touch versions.yml
"""
stub:
"""
touch versions.yml ${meta.id}.pdf
"""
}
