scatterplot <- function(out.path=getwd(), core.name, name,
    intensity, MAnorm, array1.name, array2.name, cex, graph,
    marray, fast=50000) {

  # Output file name.
  file.out <- .escape(
      paste(core.name,'_XY_', .dtag(), '.png',sep=''));
  file.out <- file.path(out.path, file.out);

  # Keep only mapped probes.
  MAnorm <- MAnorm[MAnorm$probeID %in% marray$mapping$probeID,];
  if (fast) {
    MAnorm <- MAnorm[sample(1:nrow(MAnorm), size=fast),];
  }

  xlim <- ylim <- range(
      MAnorm[,grep('^M[12].norm$', colnames(MAnorm))]
  );

  if (graph == 'x11') {
    png(file=file.out, width=400, height=450);
    cex.point=cex;
  }
  if (graph == 'ps')  {
    bitmap(file=file.out, width=400, height=450, point=1,
        units='px', taa=4);
    cex.point=0.01*cex;
  }

  par(cex=cex);

  # Make the main plot name.
  main <- paste(name, sub('([^%])$', '\\1%', intensity));
  # Plot an empty frame.
  plot(
      xlim,
      ylim,
      type = 'n',
      main = main,
      cex.main = cex,
      cex.lab = 0.7 * cex,
      xlab = array1.name,
      ylab = array2.name,
      pch = '.',
      cex = cex
  );

  # Do the scatter plot.
  points(MAnorm$M1.norm, MAnorm$M2.norm, cex=cex.point);

  # Print version control
  .print.vcontrol("left", cex=cex);

  # Print correlation in the top-right corner.
  r <- cor(MAnorm$M1.norm, MAnorm$M2.norm, use='pairwise.complete.obs');
  text(x=xlim[2], y=ylim[2], labels=paste("r =", round(r,2)),
      adj=1, cex=cex);

  # Good bye :-)
  dev.off();

}
