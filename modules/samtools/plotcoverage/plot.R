
T2 <- read.table("TMP/depths.tsv", header = TRUE, sep = "\t", colClasses = c("character","character","numeric"),stringsAsFactors=FALSE)
head(T2)
maxDP <- max(T2$maxDP)
samplesRGB <- rainbow(nrow(T2),alpha=0.3)

# read exon bounds
T3 <- read.table("TMP/exons.bed", header = TRUE, sep = "\t", colClasses = c("character","integer","integer"),stringsAsFactors=FALSE)
head(T3)

pdf(paste("__CONTIG__","_","__START__","_","__END__",".pdf",sep=""),width=14)

# empty graphics
plot(1, type="l", main="Coverage __CONTIG__:__START__-__END__",sub="${meta.prefix?:""}",xlab="Genomic position on __CONTIG__", ylab="DEPTH", xlim=c(__START__, __END__), ylim=c(0,maxDP))

# plot exons
exonRGB <- rgb(0,0,1.0,alpha=0.05)
for (i in 1:nrow(T3)) {
 rect(xleft=T3[i,]$start,xright=T3[i,]$end,ybottom=0,ytop=maxDP,col=exonRGB,border=NA)
 }

# plot each coverage
for (i in 1:nrow(T2)) {
	T1<- read.table(T2[i,]$DEPTH, header = FALSE, sep = "\t",col.names=c('POS','DP'), colClasses = c("integer","numeric"),stringsAsFactors=FALSE)
	T1$DP <- runmed(T1$DP,k=21)
	lines(x=T1$POS,y=T1$DP,col=samplesRGB[i])
}

legend("bottomright",legend=T2$SN,title="Samples",pch=16,col=samplesRGB) 
dev.off()