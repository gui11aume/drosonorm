NGprocess <- function(
      files.table,
      out.path = paste('out', .dtag(), sep="_"),
      platform = "GPL8471",
      release = c("dm3R5", "dm2R4"),
      GFF = TRUE,
      plotMA = TRUE,
      plotXY = TRUE,
      plotACF = TRUE,
      plotProfile = TRUE,
      graph = "ps",
      cex = ifelse(graph == "ps", 1.8, 1)
) {
# Important: assume array names are of the form "CMF_123456_..."

  # Array specifications. Try to lazy-load data.
  release <- match.arg(release);
  array.object <- paste(platform, release, sep="."); # Ex: GPL8471.dm3R5
  if (exists(array.object)) {
    assign("marray", get(array.object));
  }
  else {
    stop(paste("array", array.object, "does not exist"));
  }

  # Read-in the pair table in 'meta' data.frame
  meta <- read.delim(files.table, colClasses="character",
      comment.char="#", stringsAsFactors=FALSE);

  # Make sure all required columns are present in the meta file
  # and turn them to lowercase.
  colnames(meta) <- .checknames(colnames(meta)) # Convenience function.


  # Check that the pair files exist.
  # Look for column names of the form "Exp1", "Ctl2" etc.
  Exp.ind <- grep("^exp[0-9]+$", colnames(meta));
  Ctl.ind <- grep("^ctl[0-9]+$", colnames(meta));
  allpairfiles <- c(unlist(meta[,Exp.ind]), unlist(meta[,Ctl.ind]));
  notfound <- which(allpairfiles != "" & ! file.exists(allpairfiles));
  # Stop if any pair file is not found.
  if (any(notfound)) {
    if (length(notfound) > 20) {
      # If the string passed to stop() is too long, it is truncated.
      # Truncate manually and add "..." if more than 20 files not found.
      notfound <- c(notfound[1:20], length(allpairfiles)+1);
      allpairfiles <- c(allpairfiles, "...");
    }
    string <- paste(allpairfiles[notfound], collapse="\n   ");
    stop (paste(sep="\n   ", "The following file(s) were not found:",
        paste(allpairfiles[notfound], collapse="\n   ")));
  }


  # Create out.path if it does not exist.
  if (!file.exists(out.path)) {
    dir.create(out.path);
    base::cat(paste("\nwriting output files in", out.path, "\n"));
  }

  # Main loop: run over the lines of meta.
  for (i in 1:nrow(meta)) {

    # Get the protein name.
    protname <- meta$name[i];

    # Count the arrays by the presence of "Exp2" data.
    n.arrays <- ifelse(meta$exp2[i] == "", 1, 2);

    # Get array names.
  	 array1.name <- sub("CMF_([0-9]+)_.*","\\1", meta$exp1[i]);

  	 if (n.arrays == 2) {
      array2.name <- sub("CMF_([0-9]+)_.*","\\1", meta$exp2[i]);
    }
    else {
      array2.name <- ""; # Empty string if not present.
    }

    # Make "core.name" for file naming.
    # Example core name for two arrays: HP1_70.
    core.name <- paste(
        protname,
        gsub("%", "", meta$intensity[i]), # Remove '%' sign.
        sep = "_"
    );
    # Add array name if only one.
    if (n.arrays == 1) {
      # Example core name for one array: HP1_331836_70.
      core.name <- sub(
          "_",
          paste("_", array1.name, "_", sep = ""),
          core.name
      );
    }

    base::cat(paste("\n-- processing:", meta$name[i],
        meta$intensity[i], "\n\n"));

    ### Normalization (always).
    base::cat("normalizing...\n");
    MAnorm <- NULL;
    try(expr = 
       MAnorm <- loess.norm(
           exp1.file = meta$exp1[i],
           ctl1.file = meta$ctl1[i],
           exp2.file = meta$exp2[i],
           ctl2.file = meta$ctl2[i],
           marray = marray
      )
    );
    # Write norm data.frame to file.
    if (!is.null(MAnorm)) {
      # Output file name.
      out.file <- .escape(
          paste(core.name, "_norm_", .dtag(), ".txt", sep=""));
      out.file <- file.path(out.path, out.file);
      write.table(
          MAnorm,
          file = out.file,
          row.names = FALSE,
          quote = FALSE,
          sep = "\t"
      );
    }

    # Format to gff.
    if (GFF) {
      base::cat("creating gff file...\n");
      gff <- NULL;
      try(expr =
         gff <- norm2gff(
             name = protname,
             core.name = core.name,
             MAnorm = MAnorm,  # Just created.
             marray = marray
         )
      );
      # Write gff data.frame to file.
      if (!is.null(gff)) {
        # Output file name.
        out.file <- .escape(
             paste(core.name, "_", .dtag(), ".gff", sep=""));
        out.file <- file.path(out.path, out.file);
        write.table(
            gff,
            file = out.file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
        );
      }

    }

    ### various plots
    if (plotMA) {
      base::cat("creating MA plot...\n");
      try(expr =
         plot.MA(
             out.path = out.path,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             norm = MAnorm,
             array1.name = array1.name,
             array2.name = array2.name,
             exp1spike = meta$exp1spike[i],
             ctl1spike = meta$ctl1spike[i],
             exp2spike = meta$exp2spike[i],
             ctl2spike = meta$ctl2spike[i],
             marray = marray,
             cex = cex,
             graph = graph
         )
      );
    }

    if (plotXY && (n.arrays == 2)) {
      base::cat("creating scatter plot...\n");
      try(
         scatterplot(
             out.path = out.path,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             norm = MAnorm,
             array1.name = array1.name,
             array2.name = array2.name,
             cex = cex,
             graph = graph,
             marray = marray
         )
      );
    }

    if (plotACF) {
      base::cat("creating acf plot...\n");
      try(expr =
         plot.acf.norm(
             out.path = out.path,
             marray = marray,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             norm = MAnorm,
             array1.name = array1.name,
             array2.name = array2.name,
             maxlag = 200,
             cex = cex,
             graph = graph
        )
      );
    }

    if (plotProfile) {
      base::cat("creating profiles...\n");
      try(expr =
         plot.profile(
             out.path = out.path,
             core.name = core.name,
             name = protname,
             gff = gff,
             intensity = meta$intensity[i],
             cex = cex,
             graph = graph
        )
      );
    }
  }
}
