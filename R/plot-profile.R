plot.profile <- function(core.name, name, gff,
    intensity, cex, graph) {

  out.file <- paste(core.name, '_profiles_', .dtag(), '.png', sep='');

  if (graph == 'x11') {
    png(file=out.file, width=1000, height=600);
    cex.graph=cex;
    par(mar=c(5,4,2,2));
  }
  if (graph == 'ps') {
    bitmap(file=out.file, width=1000, height=600, units='px', point=1,
        taa=4);
    cex.graph = 1.5*cex;
    par(mar=c(6,5,2,2));
  }

  par(mfrow=c(2,2), bty="n");


  ##########################################
  ##         'DRY' plot function          ##
  ##########################################

  # Important for 'plotpanel()' to work.
  colnames(gff) <- c("seqname", "source", "feature", "start", "end",
      "score", "strand", "frame", "attribute");
  ylim <- quantile(gff$score, probs=c(0.001,0.999), na.rm=TRUE);

  # Make space for genes below the profile.
  ylim[1] <- ylim[1] - 0.2*(ylim[2]-ylim[1]); 
  ysize <- ylim[2]-ylim[1];

  plot.panel <- function(chrom, xlim) {
    # Get probes on given chromosome, and then in given region.
    values <- gff[gff$seqname == chrom,];
    values <- gff[gff$start > xlim[1] & gff$end < xlim[2],];

    # Get the genes of given chromosome ('genesdef' is lazy-loaded),
    # then get genes in the given region, then separate plus and minus.
    genes <- genesdef[genesdef$seqname == chrom,];
    genes <- genes[genes$start > xlim[1] & genes$end < xlim[2],];
    genes.plus <- genes[genes$strand == "+",];
    genes.minus <- genes[genes$strand == "-",];

    # Plot (x axis is in kb).
    plot(
        x = values$start/1000,
        y = values$score,
        type = "h",
        cex.axis = 0.5 * cex.graph,
        cex.lab = cex.graph,
        xlab = paste(chrom, "- position (kb)"),
        ylab = "log2(Dam fusion/Dam)",
        ylim = ylim,
        lwd = 1
    );

    # Plot the genes.
    rect(
        xleft = genes.plus$start/1000,
        xright = genes.plus$end/1000,
        ybottom = ylim[1] + 0.05 *ysize,
        ytop = ylim[1] + 0.1 * ysize,
        col = "white",
        border = TRUE,
        lwd=0.2
    );
    rect(
        xleft = genes.minus$start/1000,
        xright = genes.minus$end/1000,
        ybottom = ylim[1],
        ytop = ylim[1] + 0.05 * ysize,
        col = "black",
        border = TRUE,
        lwd = 0.2
    );
  }

  xlim <- c(1, 2.5)*1e6; # pericentric 2R
  plot.panel(chrom='chr2R', xlim=xlim); 
  # Write the protein name in the top-left corner.
  text(x=xlim[1]/1000, y=ylim[2], paste(name, intensity), col='grey50',
      cex=cex.graph*2, adj=c(0,0.5));

  plot.panel(chrom='chr2L', xlim=c(13.75, 13.95)*1e6); # 2L hi-density

  # Print version control.
  .print.vcontrol("right", cex=cex);

  plot.panel(chrom='chr4', xlim=c(0, 1.281)*1e6);      # whole chr4 
  plot.panel(chrom='chr3R', xlim=c(12.35, 13)*1e6);    # BX-C

  dev.off();

}
