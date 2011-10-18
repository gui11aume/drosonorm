### Loess NG arrays and make MA plots, scatterplots and autocorrelation plots
# Editor: Guillaume Filion (modifications to build a pipeline on Zircon)
# Global edits: 'T' and 'F' converted to 'TRUE' and 'FALSE'.
# NOTE: I prepended 'Guillaume Filion:' to the edits.
# Date: June 28, 2011.
# U. Braunschweig, 10/10/2009
# Functions for loess normalization, GFF conversion and various plots
# to be called by 'NGprocess()'
#
# Arrays are separately loess normalized, then median centered and averaged.
# GFF comprises all probes on 'chromosomes', except those with unspecified start
# (this is for compatibility with liftover annotations) and those overlapping with
# fusion protein ORF (if desired and DID specified)
# 
# Changes: zeroes in pair files that would cause Inf or -Inf M/A values are set to NA; NAs handled in plots

# Guillaume Filion: functions in the script:
# single.loess()
# NormToGFF()
# Draw.MA.norm()
# Draw.MA.single()
# spikes2bin()
# Draw.XY.norm()
# Draw.ACF.norm()
# Draw.ACF.single()
# Draw.profiles()

#### loess normalization of 1 array
#single.loess <- function(exp1.file,ctl1.file, span=0.75, core.name, out.path){
#  cat('Doing normalization...   ')
# 	array1.name <- sub('[CMF_]*([0-9]+)_.*','\\1',as.character(exp1.file))
# 	
#	exp1 <- read.delim(exp1.file, comment.char='#',colC=c(rep('character',4),rep('numeric',7)))
#	ctl1 <- read.delim(ctl1.file, comment.char='#',colC=c(rep('character',4),rep('numeric',7)))
#	exp1 <- exp1[order(exp1$PROBE_ID),]
#	ctl1 <- ctl1[order(ctl1$PROBE_ID),]
#	setwd(out.path)
#	
#	ma <- cbind(exp1[,c(3,4,10)], ctl1[,10])
#	names(ma) <- c("SEQ_ID", "PROBE_ID", "exp1", "ctl1")
#	ma <- cbind(ma, M=log2(ma$exp1/ma$ctl1), A=0.5*log2(ma$exp1*ma$ctl1))
#  fit <- rep(NA,nrow(ma))
#	fit[!is.na(ma$M)] <- loess(M ~ A, ma, span=span, control=loess.control(trace.hat='approximate'), na.action=na.omit)$fitted	# Do normalization
#	norm <- cbind(ma[,1:2], M=ma$M-fit)
#  norm$M <- round(norm$M-median(norm$M, na.rm=TRUE),digits=3)
#  cat('done.\n')
#
#	write.table(ma[,2:6], file=paste(core.name, "_raw.txt", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
#	write.table(norm[,2:3], file=paste(core.name, "_norm.txt", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
#}




### MA plots for 1 single array, raw and normalized
#Draw.MA.single <- function(core.name, DID, name, intensity, raw, norm, ann, array1.name, Exp1spike, Ctl1spike, cex, graph){
#  cat('Drawing MA plots...      ')
#
#	c.between <- paste("Array 1 raw - norm.: ", round(cor(raw$M, norm$M, use='p'), digits=2), sep="")
#	c. <- matrix(nrow=6, ncol=1)	#lookup control, empty and data features
#
#	for (i in (1:6)) {  # vector with IDs of the 6 spike-in oligos
#		c.[i] <- max(0,which(ann$PROBE_ID == paste("DMEL000041096", i-1, sep="")))
#	}
#	rand. <- which(ann$SEQ_ID == "RANDOM")
#	maxx <- max(c(max(raw$A, norm$A)), na.rm=TRUE)
#	minx <- min(c(min(raw$A, norm$A)), na.rm=TRUE)
#	maxy <- max(c(max(raw$M, norm$M)), na.rm=TRUE)
#	miny <- min(c(min(raw$M, norm$M)), na.rm=TRUE)
#
#  # Test spike presence and prepare strings for plotting  
#  spikes1 <- spikes2bin(Exp1spike)-spikes2bin(Ctl1spike)
#  spikes1.hi  <- paste(sub('1','o',sub('0',' ',sub('-1',' ',spikes1))),sep='',collapse='')
#  spikes1.mid <- paste(sub('1',' ',sub('0','o',sub('-1',' ',spikes1))),sep='',collapse='')
#  spikes1.lo  <- paste(sub('1',' ',sub('0',' ',sub('-1','o',spikes1))),sep='',collapse='')
#
#  # Draw MA plots to file; random oligos are blue, the 6 spiking oligos are red
#  if (graph == 'x11') {
#    png(file=paste(core.name,'_MA.png',sep=''), wid=341, hei=768)
#    cex.point=cex
#  }
#  if (graph == 'ps')  {
#    bitmap(file=paste(core.name,'_MA.png',sep=''), wid=341, hei=768, units='px', point=1, taa=4, fonts=c('sans','mono'))
#    cex.point=0.01*cex
#  }
#  	par(mfrow=c(2,1), mar=c(5,4,3,1), pch=".", cex=cex)
#  	ylim <- c(miny,maxy)
#  	xlim <- c(minx,maxx)
#  	plot(1, 1, type='n', xlim=xlim, ylim=ylim, xlab='A', ylab='M', main=paste("Array", array1.name, intensity))  # Array 1 raw
#      points(raw$A, raw$M, cex=cex.point) 
#    	points(raw$A[rand.], raw$M[rand.], pch=20, col="blue")
#    	points(raw$A[c.], raw$M[c.], pch=20, col="red")
#    	for (i in (1:6)) {
# 	      if (max(c.) == 0) {break}
#    		text(raw$A[c.[i]]+0.1, raw$M[c.[i]], i, adj=0, col="red")
#    	}
#    	points(raw$A, raw$M-norm$M, col="red")   
# 	    par(family='mono')
#      legend(x='topright',legend=c('Spikes: 123456',paste(c('    hi: ','   mid: ','    lo: '),c(spikes1.hi,spikes1.mid,spikes1.lo),sep="")), cex=0.9, text.col='grey50',bty='n')
#      par(family='sans')
#  	plot(1, 1, type='n', xlim=xlim, ylim=ylim, xlab='A', ylab='M',, main=paste('Array', array1.name, 'normalized'))  # Mean 1+2 normalized
#      points(raw$A, norm$M, cex=cex.point)  # Mean 1+2 normalized
#    	points(raw$A[rand.], norm$M[rand.], pch=20, col="blue")
#    	text(xlim[2], maxy, c.between, adj=1)
#      legend(x='topleft',legend=c(name,DID), cex=1.5, bty='n', text.col='grey50')
# 	dev.off() 
#  cat('done.\n')                              
#}


#### Draw autocorrelation plot for 1 single array
#Draw.ACF.single <- function(core.name, DID, name, intensity, norm, array1.name, cex, graph){
#  cat('Drawing ACF plot...      ')
#  lag1 <- 200
#  ylim=c(0,1)
#  if (graph == 'x11') {png(paste(core.name,'_ACF.png',sep=''), wid=300, hei=400)}
#  if (graph == 'ps')  {bitmap(paste(core.name,'_ACF.png',sep=''), wid=300, hei=400, units='px', point=1, taa=4)}
#    par(mfrow=c(1,1),cex=cex)
#    acf.1 <- acf(norm$M,lag.max=lag1,na.action=na.pass,ci=0,main=paste('ACF of array',array1.name),ylim=ylim,col='red')
#      lines(c(-20,250),c(0.2,0.2),col='grey',lty=2)
#      lines(c(-20,250),c(0.4,0.4),col='grey',lty=2)
#      lines(c(-20,250),c(0.6,0.6),col='grey',lty=2)
#      legend(x='topright',legend=c(DID,name,intensity), cex=2, bty='n')
#  dev.off()
#  cat('done.\n')
#}


