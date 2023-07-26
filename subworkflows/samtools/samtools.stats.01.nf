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

include {SAMTOOLS_SAMPLES02} from '../../subworkflows/samtools/samtools.samples.02.nf'
include {moduleLoad;isBlank} from '../../modules/utils/functions.nf'
include {SAMTOOLS_STATS_01 as ST_STATS} from '../../modules/samtools/samtools.stats.01.nf'
include {COLLECT_TO_FILE_01} from '../../modules/utils/collect2file.01.nf'
include {MULTIQC_01} from '../../modules/multiqc/multiqc.01.nf' 
include {MERGE_VERSION} from '../../modules/version/version.merge.02.nf'
include {SIMPLE_ZIP_01} from '../../modules/utils/zip.simple.01.nf' addParams(compression_level:9)

workflow SAMTOOLS_STATS_01 {
	take:
		meta
		bams
		sample2pop /* or NO_FILE */
	main:
		version_ch = Channel.empty()
		
		samples_ch = SAMTOOLS_SAMPLES02( [with_header:true,allow_multiple_references: true,allow_duplicate_samples : true], "", bams)
		version_ch = version_ch.mix(samples_ch.version)

		stats_ch = ST_STATS([gzip:false],samples_ch.output.splitCsv(header:true,sep:"\t").map{T->T.plus("sample":T.new_sample)})
		version_ch = version_ch.mix(stats_ch.version)


		st_stats_outputs_ch = stats_ch.output.map{T->T[1]}



		if(sample2pop.name.equals("NO_FILE")) {
			toqc =  st_stats_outputs_ch
			}
		else
			{
			statsperpop_ch = STATS_PER_POP([:], stats_ch.output.map{T->T[0].sample+"\t"+T[1]}.collect(), sample2pop)
			version_ch = version_ch.mix(statsperpop_ch.version)

			toqc =  st_stats_outputs_ch.mix(statsperpop_ch.output.flatten())
			}

                file_list_ch = COLLECT_TO_FILE_01([:], toqc.collect() )
                version_ch = version_ch.mix(file_list_ch.version)


                multiqc_ch = MULTIQC_01([extra:" --fullnames "],file_list_ch.output)
                version_ch = version_ch.mix(multiqc_ch.version)

		to_zip = Channel.empty().mix(st_stats_outputs_ch).mix(multiqc_ch.zip)
		
		ch1_ch = SIMPLE_ZIP_01(to_zip.collect())
		
		version_ch = MERGE_VERSION("Samtools stats",version_ch.collect())
	emit:
		version = version_ch
		multiqc_zip = multiqc_ch.zip
		files = file_list_ch.output
		zip = ch1_ch.zip
	}

process STATS_PER_POP {
tag "N=${L.size()} ${sample2pop.name}"
input:
	val(meta)
	val(L)
	path(sample2pop)
output:
	path("*_mqc.png"),emit:output
	path("version.xml"),emit:version
script:
"""
${moduleLoad("r")}
mkdir -p TMP
cat << EOF | sort -T TMP -t '\t' -k1,1 > TMP/jeter.a
${L.join("\n")}
EOF
sort -T TMP -t '\t' -k1,1 '${sample2pop}' > TMP/jeter.b

# pop and STATS
join -t '\t' -1 1 -2 1 -o '2.2,1.2' TMP/jeter.a TMP/jeter.b > TMP/jeter.c
test -s TMP/jeter.c

cat << '__EOF__' > TMP/jeter.R

pop2file <- read.table("TMP/jeter.c",sep='\t',header=FALSE,stringsAsFactors=FALSE,colClasses = c("character","character"),col.names=c("POP","STATS"))


populations <- sort(unique(pop2file\$POP))
colors <- rainbow(length(names(populations)))


plotStats <- function(meta) {
all_collections <- c()

pop2values <- list()

i <- 1
for(pop in populations ) {
write(paste("#POP=",pop), stderr())

#initialize for pop
statsfiles <- pop2file[pop2file\$POP==pop,]\$STATS
write(statsfiles,stderr())

hash <- list()
hash\$titles <- paste(pop," N=",as.character(length(statsfiles)))
hash\$color <- colors[i]

vals<-c()

for(statf in statsfiles) {
	write(paste("#POP=",pop," STATS:",statf), stderr())

	T1 <- read.delim(statf,sep='\t',header=FALSE,stringsAsFactors=FALSE)

	T1 <- T1[T1\$V1=="SN" & T1\$V2 == meta\$column ,]
	T1 <- as.numeric(T1\$V3)
	vals<-c(vals,T1)
	}
hash\$data <-  vals
all_collections = c(all_collections,list(hash))
i <- i+1
}

png(paste(meta\$basename,"_mqc.png",sep=""))
boxplot(
	lapply(all_collections,function(item){return(item\$data)} ),
	col=unlist(lapply(all_collections,function(item){return(item\$color)} )),
	names=unlist(lapply(all_collections,function(item){return(item\$titles)} )),
	main=meta\$main,
	ylab=meta\$ylab,
	las=2,
	cex.axis=0.5
	)
dev.off()
}


plotStats(list(
	column="error rate:",
	main="Error Rate",
	ylab="Percent",
	basename="error_rate"
	))

plotStats(list(
	column="insert size average:",
	main="Insert Size",
	ylab="size",
	basename="insert_size"
	))


plotStats(list(
	column="average quality:",
	main="Average quality",
	ylab="QUAL",
	basename="average_quality"
	))

__EOF__

R --vanilla < TMP/jeter.R

##################
cat << EOF > version.xml
<properties id="${task.process}">
        <entry key="name">${task.process}</entry>
        <entry key="description">plot samtools stats per pop</entry>
</properties>
EOF
"""
}
