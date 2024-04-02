
T2 <-read.table("__SAMPLESHEET__",
	header=FALSE,
	sep="\t",
	stringsAsFactors =FALSE,
	colClasses=c("character","character","character"),
	col.names=c("sample","condition","wig")
	)


head(T2)

capy=0.0001
for (i in 1:nrow(T2)) {
	T1 <-read.table(paste0("TMP/depth.",i,".txt"), header = FALSE, sep = "\t",colClasses = c("integer","numeric"),col.names=c("pos","cov"), stringsAsFactors=FALSE)
	if(nrow(T1) > 21 ) {
		T1$cov <-  runmed(T1$cov,k=21)
		}
	capy= max(capy, T1$cov  );
	}

capy


colors<-rainbow(length(unique(T2$condition)))
pigments <- colors[as.factor(T2$condition)]

pdf("TMP/jeter.pdf", 15, 5)

T1 <-read.table("TMP/depth.1.txt", header = FALSE, sep = "\t",colClasses = c("integer","numeric"),col.names=c("pos","cov"), stringsAsFactors = FALSE)
if(nrow(T1) > 21 ) {
	T1$cov <-  runmed(T1$cov,k=21)
	}


head(T1)

plot(x=T1$pos,
	 y=T1$cov,
         main="Wig __CHROM__:__START__-__END__ extended to __CHROM__:__XSTART__-__XEND__",
         xlab="Position",
         ylab="Wig Depth",
         xlim = c(__XSTART__,__XEND__),
         ylim = c(0,capy),
         col=  pigments[1],
         type="l",
         pch=19
	 )

for (i in 2:nrow(T2)) {
	T1 <-read.table(
		paste0("TMP/depth.",i,".txt"),
		header = FALSE,
		sep = "\t",
		colClasses = c("integer","numeric"),
		col.names=c("pos","cov"),
		stringsAsFactors = FALSE
		)
	if(nrow(T1) > 21 ) {
		T1$cov <-  runmed(T1$cov,k=21)
		}

	lines(T1$pos, T1$cov, col = pigments[i] );
	}

        
abline(v=__START__, col="green")
abline(v=__END__, col="green")

legend("topleft",
	col= pigments,
	legend=paste0(T2$sample," ",T2$condition),
	cex=0.5 ,
	text.col= pigments,
	lty=1
	)

dev.off()
