### Draw QC plots of NG data
# changed from E.d.Wit 2006
# changes: can use x11 or postscript device

falseColPlots<-function(cy3.file,cy5.file,graph='x11'){
	read.delim(cy3.file, comm='#')->cy3
	read.delim(cy5.file, comm='#')->cy5
	
	cy3.name<-sub("\\..+","",as.character(cy3.file))
	cy5.name<-sub("\\..+","",as.character(cy5.file))
	
	#make raw intensity plots
	makeIntensityPlot(cy3,cy3.name,graph)
	makeIntensityPlot(cy5,cy5.name,graph)

}	

makeIntensityPlot<-function(data, name, graph){
	im<-matrix(nrow=768,ncol=1024) #the dimensions of the array
	for(i in 1:dim(data)[1]){
		im[data$X[i],data$Y[i]]<-data$PM[i]
	}

  if (graph == 'x11') {
    png(filename=paste(name, "_intensity_plot.png", sep=""), wid=800, hei=600)
  }
  if (graph == 'ps')  {
    bitmap(file=paste(name, "_intensity_plot.png", sep=""),type='png256', wid=800, hei=600, units='px')
  } 
	image(log2(im), col=topo.colors(1000), main=name, axes=F)
	dev.off()
	
}	

