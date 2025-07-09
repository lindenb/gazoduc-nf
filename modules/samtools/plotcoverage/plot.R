
T2 <- read.table("TMP/depths.tsv", header = TRUE, sep = "\t", colClasses = c("character","character","numeric"),stringsAsFactors=FALSE)
head(T2)
maxDP <- max(T2$maxDP)
samplesRGB <- rainbow(nrow(T2),alpha=0.3)

# read exon bounds
T3 <- read.table("TMP/exons.bed", header = TRUE, sep = "\t", colClasses = c("character","integer","integer"),stringsAsFactors=FALSE)
head(T3)

pdf(paste("__PREFIX__",".pdf",sep=""),width=14)

# empty graphics
plot(1, type="l", 
	main="__CONTIG__:__START__-__END__ __TITLE__",
	sub="extended to __CONTIG__:__XSTART__-__XEND__",
	xlab="Genomic position on __CONTIG__",
	ylab="DEPTH",
	xlim=c(__XSTART__, __XEND__),
	ylim=c(0,maxDP)
	)


if(nrow(T3) > 0) {
	# plot exons
	exonRGB <- rgb(0,0,1.0,alpha=0.05)
	for (i in 1:nrow(T3)) {
	rect(xleft=T3[i,]$start,xright=T3[i,]$end,ybottom=0,ytop=maxDP,col=exonRGB,border=NA)
	}
}


runmed_k <- 1 + 2 * min((__XLEN__-1)%/% 2, ceiling(0.01*__XLEN__))
if( runmed_k < 7) {
	runmed_k <- 7
}

# plot each coverage
for (i in 1:nrow(T2)) {
	T1<- read.table(
		T2[i,]$DEPTH,
		header = FALSE,
		sep = "\t",col.names=c('POS','DP'),
		colClasses = c("integer","numeric"),
		stringsAsFactors=FALSE
		)
	T1$DP <- runmed(T1$DP,runmed_k)
	lines(x=T1$POS, y=T1$DP, col=samplesRGB[i])
}

abline(v=__START__, col="green")
abline(v=__END__, col="green")


legend("bottomright",legend=T2$SN,title="Samples",pch=16,col=samplesRGB) 
dev.off()