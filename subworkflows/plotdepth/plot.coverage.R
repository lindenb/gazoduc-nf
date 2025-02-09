	start_date <- Sys.time()
	DEBUG <- function(...) {
		cat("[DEBUG][",file=stderr())
		cat((Sys.time() - start_date), file=stderr())
		cat("]:",file=stderr())
		cat(sprintf(...), sep='', file=stderr())
		cat(".\n",file=stderr())
		}
	normalize_median <- _MEDIAN_

	DEBUG("screen : _START_ to _END_  . SV _DELSTART_ to _DELEND_");
	LARGE_SV <- _LARGE_SV_
	distance <- 1 + (_END_ - _START_)
	DEBUG("distance=%d",distance)

	smoothvalue <- 1 + 2 * min((distance-1)%/% 2, ceiling(0.01*distance))
	DEBUG("read table")
	T1 <-read.table("TMP/depth.txt", header = FALSE, sep = "\t",colClasses = c("integer","numeric"),col.names=c("pos","cov"))
	

	## smooth yaxis using runmed
	if(smoothvalue> 3 ) {
		DEBUG("runmed %d",smoothvalue)
		T1$cov <- runmed(T1$cov,smoothvalue)
		DEBUG("end runmed")
		}


	# median in outer region
	DEBUG("get median")
	# outer median
	DEBUG("N = %d",nrow(T1))
	DEBUG("N-outer = %d",nrow(T1[T1$pos < _DELSTART_ | T1$pos > _DELEND_,]))
	DEBUG("N-inner = %d",nrow(T1[T1$pos >= _DELSTART_ & T1$pos <= _DELEND_,]))
	mediany <- median(T1[T1$pos < _DELSTART_ | T1$pos > _DELEND_,]$cov)
	DEBUG("outer-mediany=%f",mediany)


	
	if( normalize_median == TRUE ) {
		DEBUG("normalize median: TRUE")
		capy <- 3.0
		if( mediany > 1 ) {
			T1$cov <- T1$cov / mediany
			# recompute median y for capy
			mediany <- median(T1[T1$pos < _DELSTART_ | T1$pos > _DELEND_,]$cov)
			DEBUG("outer-mediany=%f",mediany)
			}
	} else {
		DEBUG("normalize median: FALSE")
		}
	# inner median
	mediany_in <- median(T1[T1$pos >= _DELSTART_ & T1$pos <= _DELEND_,]$cov)
	DEBUG("inner-mediany=%f",mediany_in)

	if(normalize_median == TRUE ) {
		capy <- mediany *3.0
		} else {
		capy <- max( mediany_in, ifelse(mediany > 1.0 , mediany *3.0, 3.0))
		}
	DEBUG("capy=%f vs median=%f", capy, mediany)

	yaxis <- T1$cov

	DEBUG("plot")
	pdf("TMP/jeter.pdf", 15, 5)
	plot(x=T1$pos,y=yaxis,
	   main=paste("_SAMPLE_",ifelse(normalize_median == TRUE ," Normalized","")," Coverage _CHROM_:",format(_DELSTART_,big.mark=","),"-",format(_DELEND_,big.mark=",")," len=",format((_DELEND_-_DELSTART_),big.mark=",")," ","_TITLE_",sep=""),
	  sub="_BAM_",
	  xlab="Position",
          ylab=paste(ifelse(normalize_median == TRUE ," Normalized ",""),"Depth"),
	  xlim = c(_START_,_END_),
	  ylim = c(0,capy),
	  col=rgb(0.2,0.1,0.5,0.9) , 
	  type="l", 
	  pch=19)
	
	mediany <- ifelse(normalize_median == TRUE ,1.0, mediany )

	# outer median plot
	segments(x0=_DELSTART_,x1=_DELEND_,y0=mediany_in,col="blue")
	segments(x0=_START_,x1=_DELSTART_,y0=mediany,col="blue")
	# inner median plot
	segments(x0=_DELEND_,x1=_END_,y0=mediany,col="blue")
	abline(h=mediany*0.5, col="orange")
	abline(h=mediany*1.5, col="orange")
		

	abline(v=_DELSTART_, col="green")
	abline(v=_DELEND_, col="green")
	
	DEBUG("plot exons")
	source("TMP/exons.R",local=TRUE)
	dev.off()
	DEBUG("done")
