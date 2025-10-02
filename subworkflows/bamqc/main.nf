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
include {MOSDEPTH                   } from '../../modules/mosdepth'
include {SAMTOOLS_STATS             } from '../../modules/samtools/stats'
include {SAMTOOLS_FLAGSTATS         } from '../../modules/samtools/flagstats'
include {SAMTOOLS_IDXSTATS          } from '../../modules/samtools/bamidxstats'
include {BED_TO_INTERVAL_LIST       } from '../../modules/gatk/bed2intervallist'
include {COLLECT_MULTIPLE_METRICS   } from '../../modules/gatk/collectmultiplemetrics'

workflow BAM_QC {
take:
	meta
	fasta
	fai
	dict
	bed /* meta, bed */
	bams /* meta, bam, bai */
main:
	versions_ch = Channel.empty()
	reports_ch  = Channel.empty()

  /***************************************************
   *
   *  MOSDEPTH
   *
   */
	MOSDEPTH(
		fasta,
		fai,
		bams.combine(bed).map{[it[0],it[1],it[2],it[4]]} ,
		)
	versions_ch = versions_ch.mix(MOSDEPTH.out.versions)

	mosdepth_global  = MOSDEPTH.out.global_txt
	mosdepth_summary = MOSDEPTH.out.summary_txt
	mosdepth_regions = MOSDEPTH.out.regions_txt

	reports_ch = reports_ch
		.mix(mosdepth_global)
		.mix(mosdepth_summary)
		.mix(mosdepth_regions)

	/** PLOT chrX / chrY */
	ch1 = mosdepth_summary.splitCsv(header:true,sep:'\t')
		.filter{it[1].chrom.matches("(chr)?[XY]_region")}
		.map{[it[1].chrom.replaceAll("_region",""), it[0].id, (it[1].mean?:it[1].median)]}
		.branch{v->
			chrX: v[0].contains("X")
			chrY: v[0].contains("Y")
			}
	
	PLOT_CHR_XY(
		ch1.chrX.join(ch1.chrY, by:1) 
			.map{[it[0],it[2],it[4]]}//sample,depthX,depthY
			.map{it.join("\t")}
			.toSortedList()
			.map{[[id:"xy"],it]}
		)
	versions_ch = versions_ch.mix(PLOT_CHR_XY.out.versions)

	/** outliers of DEPTH */
	treshold_ch = Channel.of(100,1000);

	ch1 = mosdepth_summary
		.combine(treshold_ch)
		.map{meta,summary,fold->[meta.id,meta.plus(treshold:fold),summary]}
		.join(bams.map{[it[0].id,it[0],it[1],it[2]]})
		.combine(bed)
		.map{sample,meta1,summary,meta2,bam,bai,meta3,bed->[meta1,summary,bam,bai,bed]}

	DEPTH_OUTLIER(fasta,fai,ch1)
	versions_ch = versions_ch.mix(DEPTH_OUTLIER.out.versions)

	DEPTH_OUTLIER_MERGE_ALL_SAMPLES(DEPTH_OUTLIER.out.bed.map{meta,bed,tbi->[[id:"outliers",treshold:meta.treshold],bed]}.groupTuple())
	versions_ch = versions_ch.mix(DEPTH_OUTLIER_MERGE_ALL_SAMPLES.out.versions)
	
  /***************************************************
   *
   *  SAMTOOLS STATS
   *
   */
    SAMTOOLS_STATS(fasta,fai,bed,bams)
	versions_ch = versions_ch.mix(SAMTOOLS_STATS.out.versions)
	reports_ch = reports_ch.mix(SAMTOOLS_STATS.out.stats)

  /***************************************************
   *
   *  SAMTOOLS IDXSTATS
   *
   */
    SAMTOOLS_IDXSTATS(bams)
	versions_ch = versions_ch.mix(SAMTOOLS_IDXSTATS.out.versions)
	reports_ch = reports_ch.mix(SAMTOOLS_IDXSTATS.out.stats)

  /***************************************************
   *
   *  SAMTOOLS IDXSTATS
   *
   */
    SAMTOOLS_FLAGSTATS(
		fasta,
		fai,
		bed,
		bams
		)
	versions_ch = versions_ch.mix(SAMTOOLS_FLAGSTATS.out.versions)
	reports_ch = reports_ch.mix(SAMTOOLS_FLAGSTATS.out.stats)



  /***************************************************
   *
   *  GATK QC
   *
   */
	BED_TO_INTERVAL_LIST(
		dict,
		bed
		)
	versions_ch = versions_ch.mix(BED_TO_INTERVAL_LIST.out.versions)

	COLLECT_MULTIPLE_METRICS(
		fasta,
		fai,
		dict,
		[[:],[]],//refflat
		BED_TO_INTERVAL_LIST.out.interval_list,
		bams
		)
	versions_ch = versions_ch.mix(COLLECT_MULTIPLE_METRICS.out.versions)


	reports_ch = reports_ch
		.mix(COLLECT_MULTIPLE_METRICS.out.alignment_summary_metrics)
		.mix(COLLECT_MULTIPLE_METRICS.out.base_distribution_by_cycle_metrics)
		.mix(COLLECT_MULTIPLE_METRICS.out.insert_size_metrics)
		.mix(COLLECT_MULTIPLE_METRICS.out.quality_by_cycle_metrics)
		.mix(COLLECT_MULTIPLE_METRICS.out.quality_distribution_metrics)

emit:
	versions = versions_ch
	multiqc = reports_ch
	mosdepth_global
	mosdepth_summary
	mosdepth_regions
}


process PLOT_CHR_XY {
label "process_single"
tag "${meta.id?:""}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
   tuple val(meta),val(data)
output:
    tuple val(meta),path("*.tsv"),emit:tsv
	tuple val(meta),path("*.png"),emit:png
	path("versions.yml"),emit:versions
script:
    def prefix = task.ext.prefix?:"${meta.id?:"XY"}"
	def treshold  = task.ext.treshold?:10.0
	def treshold2 = task.ext.treshold2?:2.0
	def subtitle  = task.ext.subtitle?:"Treshold M/F : ${treshold}"
"""
mkdir -p TMP
cat << EOF > TMP/jeter.tsv
${data.join("\n")}
EOF

echo 'sample\tdpx\tdpy\tratio\tsex' > TMP/samples.sex.tsv


awk -F '\t' '{S="male";FX=(\$2*1.0);FY=(\$3*1.0); if(FX > (FY * ${treshold2} )) {S="ambigous"}; if(FX > (FY * ${treshold} )) {S="female"};printf("%s\t%s\t%s\\n",\$0,(FX<=0.0?".":(FY/FX)),S);}' TMP/jeter.tsv >> TMP/samples.sex.tsv


cat << '__EOF__' | R --vanilla
T1<-read.table("TMP/samples.sex.tsv",header = TRUE,sep="\t",comment.char="",stringsAsFactors=FALSE)
male <-T1[T1\$sex=="male",]
head(male)

female <-T1[T1\$sex=="female",]
head(female)

ambigous <-T1[T1\$sex=="ambigous",]
head(ambigous)


png("TMP/${prefix}.png")
plot(1,
        xlab="Depth chrX",
        ylab="Depth chrY",
        main="Sex guessed from Depth of Coverage.",
        sub="${subtitle}",
        xlim=c(0,max(T1\$dpx)),
        ylim=c(0,max(T1\$dpy))
        )


text(x = ambigous\$dpx, y = ambigous\$dpy, labels = ambigous\$sample , pos = 4 , cex = 0.5 , col= "green") 

mc <- rgb(0,0,1.0,alpha=0.5)
points(x=male\$dpx,y=male\$dpy,type='p',col=mc,pch=16,
        xlim=c(0,max(T1\$dpx)),
        ylim=c(0,max(T1\$dpy))
        )

fc <- rgb(1.0,0,0,alpha=0.5)
points(x=female\$dpx,y=female\$dpy,type='p',col=fc,pch=16,
        xlim=c(0,max(T1\$dpx)),
        ylim=c(0,max(T1\$dpy))
        )

ac <- rgb(0.5,0.8,0.5,alpha=0.5)
points(x=ambigous\$dpx,y=ambigous\$dpy,type='p',col=ac,pch=16,
        xlim=c(0,max(T1\$dpx)),
        ylim=c(0,max(T1\$dpy))
        )


legend("topright",legend=c("male","female","ambigous"),title="Sex",pch=16,col=c(mc,fc,ac)) 
dev.off()
__EOF__

mv TMP/samples.sex.tsv ./${prefix}.tsv
mv TMP/${prefix}.png ./

touch versions.yml
"""
}


process DEPTH_OUTLIER {
label "process_single"
tag "${meta.id?:""} ${meta.treshold}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
   tuple val(meta1),path(fasta)
   tuple val(meta2),path(fai)
   tuple val(meta),path(summary),path(bam),path(bai),path(optional_bed)
output:
	tuple val(meta),path("*.bed.gz"),path("*.tbi"),optional:true,emit:bed
	path("versions.yml"),emit:versions
script:
	def treshold = (meta.treshold as double)
	def prefix = task.ext.prefix?:"${meta.id}_t${treshold}"
"""
mkdir -p TMP

awk -F '\t' '(\$1!="total_region" && \$1 ~ /_region\$/) ' "${summary}" |\\
	cut -f1,4 |\\
	sed 's/_region//' | while read CONTIG DEPTH
	do

		if ${optional_bed?true:false}
		then
			awk -F '\t' -vC="\${CONTIG}" '\$1==C' "${optional_bed}" > TMP/capture.bed
			if ! test -s TMP/capture.bed
			then
				echo "\${CONTIG}\t0\t1" > TMP/capture.bed
			fi
		fi

		samtools view --uncompressed -F 3844  -T "${fasta}" "${bam}" "\${CONTIG}" |\\
			samtools depth  ${optional_bed?"-b TMP/capture.bed":""} '-' |\\
			awk -F '\t' -vDP="\${DEPTH}" '(int(\$3) > ${treshold}*DP ) {printf("%s\t%d\t%s\t%s\\n",\$1,int(\$2)-1,\$2,\$3);}' |\\
			bedtools merge -c 4 -o median >> TMP/jeter.bed
	done

if test -s TMP/jeter.bed
then
	sort -T TMP -t '\t' -k1,1 -k2,2n TMP/jeter.bed | bgzip > TMP/jeter.bed.gz
	tabix -f -p bed TMP/jeter.bed.gz
	mv TMP/jeter.bed.gz ${prefix}.bed.gz
	mv TMP/jeter.bed.gz.tbi ${prefix}.bed.gz.tbi
fi

touch versions.yml
"""
}



process DEPTH_OUTLIER_MERGE_ALL_SAMPLES {
label "process_single"
tag "${meta.id?:""} ${meta.treshold}"
afterScript "rm -rf TMP"
conda "${moduleDir}/../../conda/bioinfo.01.yml"
input:
   tuple val(meta),path("BED/*")
output:
	tuple val(meta),path("*.bed.gz"),path("*.tbi"),emit:bed
	path("versions.yml"),emit:versions
script:
	def treshold = (meta.treshold as double)
	def prefix = task.ext.prefix?:"${meta.id}_t${treshold}"
"""
mkdir -p TMP

bedtools multiinter  -i `find BED/ -name "*.bed.gz"`| bgzip > TMP/jeter.bed.gz
tabix -f -p bed TMP/jeter.bed.gz

mv TMP/jeter.bed.gz "${prefix}.bed.gz"
mv TMP/jeter.bed.gz.tbi "${prefix}.bed.gz.tbi"


touch versions.yml
"""
}
