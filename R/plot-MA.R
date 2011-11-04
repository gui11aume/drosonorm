plot.MA <- function(out.path=getwd(), core.name, name, intensity,
    MAnorm, array1.name, array2.name, exp1spike, ctl1spike,
    exp2spike="", ctl2spike="", marray, cex, graph, fast=50000) {

  # Output file name.
  out.file <- .escape(
      paste(paste(core.name, "_MA_", .dtag(), ".png", sep="")));
  out.file <- file.path(out.path, out.file);


  # Get plot limits.
  xlim <- range(MAnorm[,grep("^A$|^A[12]", names(MAnorm))], na.rm=TRUE);
  ylim <- range(MAnorm[,grep("^M.norm$|^M[12]", names(MAnorm))], na.rm=TRUE);

  if (fast) {
     MAnorm <- rbind(
         subset(MAnorm, MAnorm$probeID %in% marray$spikes),
         subset(MAnorm, 1:nrow(MAnorm) %in% sample(1:nrow(MAnorm), fast))
     );
  }

  # Get the number of arrays.
  n.arrays <- ifelse(missing(array2.name) || array2.name == "", 1, 2);

  # Add the "%" to intensity (if needed).
  intensity <- paste(sub("%$", "", intensity), "%", sep="");


  # Set plot parameters to draw MA plots to file.
  canvas.width <- switch(
      EXPR = n.arrays,
        "1" = 350,
        "2" = 1024
  );
  if (graph == "x11") {
    png(file = out.file,
        width = canvas.width,
        height = 768);
  }
  else if (graph == "ps")  {
    bitmap(file = out.file,
           width = canvas.width,
           height = 768,
           units = "px",
           taa = 4,
           point = 1,
           fonts = c("sans", "mono"));
  }

  # Design the layout.
  if (n.arrays == 1) {
    par(mfrow=c(2,1), mar=c(5,4,3,1), pch=".", cex=cex);
  }
  else {
    par(mfrow=c(2,3), mar=c(5,4,3,1), pch=".", cex=cex);
  }



  ##########################################
  ##           spike functions            ##
  ##########################################

  spikes2bin <- function(x) {
  # Translate spike names "1a"-"31b" into oligo present/absent calls
  # Forr type-a spikes, the presence/absence is the binary representation
  # of the oligo number. For type-b spikes it is the complementary.

    # Check spike specification.
    if (!grepl("^[0-9]{1,2}(a|b)$", x)) { # 2 digits, followed by a or b.
      return (NA);
    }

    # Call to "sub" is safe because the format has been checked.
    ab <- sub(".*(a|b)$", "\\1", x); 
    x <- as.integer(sub("^([0-9]{1,2}).*", "\\1", x));

    if (x > 32) {
      return (NA);
    }
    if (ab == "b") {
      x <- 63-x;
    }

    spikes <- rep(NA_integer_, 6); 
    for (i in 1:6) {
      spikes[i] <- x %% 2;
      x <- (x - spikes[i])/2;
    }

    return (spikes); # A 0-1 vector or length 6.

  }

  spike.string <- function(ExpSpike, CtlSpike) {
    # Get spike positions (-1,0,1).
    spikes <- spikes2bin(ExpSpike)-spikes2bin(CtlSpike);
    # Make a string representation (e.g. "o  o o" for spikes 1, 4 and 6).
    lo <- paste(c(" ", "o")[1 + (spikes == -1)], collapse="");
    mid <- paste(c(" ", "o")[1 + (spikes == 0)], collapse="");
    hi <- paste(c(" ", "o")[1 + (spikes == 1)], collapse="");

    # Return an inset of the form:
    # Spikes: 123456 
    #     hi: o  o o 
    #    mid:        
    #    low:  oo o  
    return(c("Spikes: 123456", paste("    hi: ", hi, " ", sep=""),
        paste("   mid: ", mid, " ", sep=""),
        paste("   low: ", lo, " ", sep="")));

  }
 
  # Make the legend.
  spike.legend.1 <- spike.string(exp1spike, ctl1spike);
  if (n.arrays == 2) {
    spike.legend.2 <- spike.string(exp2spike, ctl2spike);
  }





  ##########################################
  ##          DRY  plot function          ##
  ##########################################

  plot.panel <- function(A.name, M.name, title, legend="",
       spikes=TRUE) {

    # Plot the empty frame.
    plot(xlim, ylim, type="n", xlab="A", ylab="M", main=title);

    # Plot the probes.
    points(MAnorm[[A.name]], MAnorm[[M.name]], pch=".");

    # Plot the random probes (in blue).
    points(MAnorm[marray$random, A.name], MAnorm[marray$random, M.name],
        pch=20, col="blue");

    # Plot the legend.
    par(family="mono");
    legend(x="topright", legend=legend, text.col="grey50", bty = "n");
    par(family="sans");

    if (spikes) {
      # Plot the pikes (in red).
      points(MAnorm[marray$spikes, A.name], MAnorm[marray$spikes, M.name],
          pch=20, col="red");
      # Add text to the spikes.
      text(MAnorm[marray$spikes, A.name]+0.1, MAnorm[marray$spikes, M.name],
          adj=0, col="red");
    }
  }

  if (n.arrays == 2) {
    MAnorm$M <- (MAnorm$M1 + MAnorm$M2) / 2;
  }

  # Start plotting.
  if (n.arrays == 1) {
    # Raw.
    plot.panel(
        A.name = "A",
        M.name = "M",
        legend = spike.legend.1,
        title = paste("Array", array1.name, intensity)
    );

    # Plot the loess line (in red).
    points(MAnorm$A, MAnorm$M-MAnorm$M.norm, col="red", pch=".");

    # Print version control on the top-left corner.
    .print.vcontrol("left", cex=cex)


    # Normalized.
    plot.panel(
        A.name = "A",
        M.name = "M.norm",
        spikes = FALSE,
        title = paste("Array", array1.name, "normalized")
    );
  }
  else { 
    # Array 1 raw.
    plot.panel(
        A.name = "A1",
        M.name = "M1",
        legend = spike.legend.1,
        title = paste("Array", array2.name, intensity)
    );

    # Plot the loess line (in red).
    points(MAnorm$A1, MAnorm$M1-MAnorm$M1.norm, col="red");

    # Print version control on the top-left corner.
    .print.vcontrol("left", cex=cex)

    # Array 2 raw.
    plot.panel(
        A.name = "A2",
        M.name = "M2",
        legend = spike.legend.2,
        title=paste("Array", array2.name, intensity)
    );

    # Plot the loess line (in red).
    points(MAnorm$A2, MAnorm$M2-MAnorm$M2.norm, col="red");

    # Mean 1+2 raw.
    plot.panel(
        A.name = "A",
        M.name = "M",
        spikes = FALSE,
        title = "Mean raw"
    );

    # Array 1 normalized.
    plot.panel(
        A.name = "A1",
        M.name = "M1.norm",
        legend = spike.legend.1,
        title = paste("Array", array1.name, "normalized")
    );

    # Array 2 normalized.
    plot.panel(
        A.name = "A2",
        M.name = "M2.norm",
        legend = spike.legend.2,
        title = paste("Array", array2.name, "normalized")
    );

    # Array 1+2 normalized.
    plot.panel(
        A.name = "A",
        M.name = "M.norm", 
        spikes = FALSE,
        title = "Mean normalized"
    );
  }

  dev.off();

}
