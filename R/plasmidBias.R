plot.plasmid.bias <- function (out.path=getwd(), core.name, name, intensity,
    norm, array1.name, array2.name, marray, gene, cex, graph) {
# Plot the signal around the ORF used in DamID.

  # Output file name.
  out.file <- .escape(
      paste(paste(core.name, '_plasmid-bias_', .dtag(), '.png', sep = '')));
  out.file <- file.path(out.path, out.file);


  # Get the number of panels to plot.
  if ("M2" %in% colnames(norm)) {
    n.arrays <- 2;
    columns <- c("M1", "M2");
    arrays <- c(array1.name, array2.name);
  }
  else {
    n.arrays <- 1;
    columns <- "M";
    arrays <- array1.name;
  }

  # Set plot parameters to draw plot(s) to file.
  canvas.width <- switch(
      EXPR = n.arrays,
        '1' = 448,
        '2' = 896
  );
  if (graph == 'x11') {
    png(file = out.file,
        width = canvas.width,
        height = 512);
  }
  else if (graph == 'ps')  {
    bitmap(file = out.file,
           width = canvas.width,
           height = 512,
           units = 'px',
           taa = 4,
           point = 1,
           fonts = c('sans', 'mono'));
  }

  par(mfrow = c(1, n.arrays), cex=cex);


  # -- Do the plot(s). -- ##

  g <- .getgene(gene=gene, release=marray$release);
  mapping <- marray$mapping;

  # Focus on the chromosome of interest.
  mapping <- subset(mapping, mapping$seqname == g$seqname);
  scores <- subset(norm, norm$seqname == g$seqname);

  for (i in 1:n.arrays) {
    plot(
       x = scores$start/1e6,
       y = scores[[columns[i]]],
       type = "h",
       xlim = c(g$start-50000, g$end+50000)/1e6,
       main = paste(
           "Plasmid bias of array",
           arrays[i],
           paste("(", sub("%$", "", intensity), "%)", sep = "")),
       xlab = paste("position on ", g$seqname, " (Mb)", sep = ""),
       ylab = "Raw DamID M value"
    );
    txt <- ifelse(
        name == gene, 
        name,
        paste(name, " (", gene, ")", sep = "")
    );
    mtext(text=txt, line=0.4, col="grey50", cex=1.5*cex);

    # Add lines at gene borders.
    abline(v = c(g$start, g$end)/1e6, col = 2, lty = 2);

    if (i == 1) {
      .print.vcontrol(pos="left", cex=cex);
    }
  }


  dev.off();

}



mask.plasmid.bias <- function (MAnorm, gene, marray) {
# Set the probes that match the plasmid used for DamID to NA.

  mapping <- marray$mapping;

  # Identify the gene.
  g <- .getgene(gene=gene, release=marray$release);
  if (is.null(g)) {
    warning (paste("no mapping information for gene", gene));
    return (NULL);
  }

  # Focus on the chromosome of the gene.
  mapping <- subset(mapping, mapping$seqname == g$seqname);

  # Mask the probes.
  masked.probes <- mapping$probeID[
      (mapping$start > g$start & mapping$start < g$end) |
      (mapping$end > g$start & mapping$end < g$end)
  ];

  MAnorm$M.norm[MAnorm$probeID %in% masked.probes] <- NA;
  return(vtag(MAnorm));

}
