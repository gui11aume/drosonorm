process <- function(
      files.table,
      out.path = paste('drosonorm', .dtag(), sep="_"),
      sortby = c("name", "filetype", "none"),
      platform = "GPL8471",
      release = c("dm3R5", "dm2R4"),
      mask = TRUE,
      GFF = TRUE,
      WIG = TRUE,
      DAM = TRUE,
      targets = TRUE,
      plotBias = TRUE,
      plotMA = TRUE,
      plotXY = TRUE,
      plotACF = TRUE,
      plotSample = TRUE,
      graph = "ps",
      cex = ifelse(graph == "ps", 1.8, 1)) {

  sortby <- match.arg(sortby);

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

  # Remove blank lines.
  meta <- meta[
      colSums(!apply(meta, 1, grepl, pattern="^[[:space:]]*$")) > 0
  ,];

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


  # Turn off scientific notation for all outputs.
  op <- options();
  on.exit(expr = options(op));
  options(scipen=10);


  # Main loop: run over the lines of meta.
  for (i in 1:nrow(meta)) {

    # Get the protein name and gene.
    protname <- meta$name[i];


    # Count the arrays by the presence of "Exp2" data.
    n.arrays <- ifelse(meta$exp2[i] == "", 1, 2);

    # Get array names.

    
    array1.name <- ifelse(
        grepl("^CMF_[0-9]+", meta$exp1[i]),
        sub("CMF_([0-9]+)_.*","\\1", meta$exp1[i]),
        paste(meta$exp1[i], meta$name[i], sep = "_")
    );

  	 if (n.arrays == 2) {
      array2.name <- ifelse(
          grepl("^CMF_[0-9]+", meta$exp2[i]),
          sub("CMF_([0-9]+)_.*","\\1", meta$exp2[i]),
          paste(meta$exp2[i], meta$name[i], sep = "_")
      );
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

    # Try to find the gene by FBgn identifier.
    FBpattern <- grep("FBgn[0-9]{7}", meta[i,], value=TRUE);
    if (length(FBpattern) == 1) {
       gene <- sub(".*(FBgn[0-9]{7}).*", "\\1", FBpattern)[1];
       base::cat("detected gene ID", gene, "\n");
       this.mask <- mask;
       this.plotBias <- plotBias;
    }
    else {
       base::cat("no gene ID: skipping plasmid bias\n");
       this.mask <- FALSE;
       this.plotBias <- FALSE;
    }

    ## -- Normalization (always). --
    base::cat("normalizing...\n");
    MAnorm <- NULL;
    try(expr = 
       MAnorm <- loess.norm(
           exp1.file = meta$exp1[i],
           ctl1.file = meta$ctl1[i],
           exp2.file = meta$exp2[i],
           ctl2.file = meta$ctl2[i],
           array1.name = array1.name,
           array2.name = array2.name,
           marray = marray,
           gene = gene,
           mask = this.mask
      )
    );
    # Write norm data.frame to file.
    if (!is.null(MAnorm)) {
      # Output file name.
      this.out.path <- switch(
          EXPR = sortby,
              "name" = file.path(out.path, protname),
              "filetype" = file.path(out.path, "norm"),
              "none" = out.path
      );
      if (!file.exists(this.out.path)) {
        dir.create(this.out.path);
      }
      out.file <- .escape(
          paste(core.name, "_norm_", .dtag(), ".txt", sep=""));
      out.file <- file.path(this.out.path, out.file);
      # Add the 'meta' line to the vheader.
      addcomment(MAnorm, "meta data", paste(meta[i,], collapse=" "));
      write.table(
          MAnorm,
          file = out.file,
          row.names = FALSE,
          quote = FALSE,
          sep = "\t"
      );
    }
    else {
      base::cat("failed to normalize: skipping\n");
      next;
    }

    # gff format.
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
        this.out.path <- switch(
            EXPR = sortby,
                "name" = file.path(out.path, protname),
                "filetype" = file.path(out.path, "gff"),
                "none" = out.path
        );
        if (!file.exists(this.out.path)) {
          dir.create(this.out.path);
        }
        out.file <- .escape(
             paste(core.name, "_", .dtag(), ".gff", sep=""));
        out.file <- file.path(this.out.path, out.file);
        # Add the 'meta' line to the vheader.
        addcomment(gff, "meta data", paste(meta[i,], collapse=" "));
        write.table(
            gff,
            file = out.file,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
        );
      }
      else {
        base::cat("  failed!\n");
      }
    }

    # wig format.
    if (WIG) {
      base::cat("creating wig file...\n");
      wigdf <- NULL;
      try(expr =
         wigdf <- norm2wig(
             MAnorm = MAnorm,  # Just created.
             marray = marray
         )
      );
      # Write wig data.frame to file.
      if (!is.null(wigdf)) {
        # Output file name.
        this.out.path <- switch(
            EXPR = sortby,
                "name" = file.path(out.path, protname),
                "filetype" = file.path(out.path, "wig"),
                "none" = out.path
        );
        if (!file.exists(this.out.path)) {
          dir.create(this.out.path);
        }
        out.file <- .escape(
             paste(core.name, "_", .dtag(), ".wig", sep=""));
        out.file <- file.path(this.out.path, out.file);
        # Add the 'meta' line to the vheader.
        addcomment(wigdf, "meta data", paste(meta[i,], collapse=" "));
        # Explicitly write vheaer.
        base::cat(vheader(wigdf), file = out.file);
        # Write the track definition line.
        base::cat(
          paste(
            "track type=wiggle_0 name=",
            meta$name[i],
            " description=DamID profile of ",
            meta$name[i],
            " (probe-wise signal, release ",
            marray$release,
            ", microarray platform ",
            marray$name, ")\n",
            sep = ""
          ),
          file = out.file,
          append = TRUE
        );
        for (seqname in unique(wigdf$seqname)) {
          base::cat(
            paste("variableStep chrom=", seqname, "\n", sep = ""),
            file = out.file,
            append = TRUE
          ); 
          write.table(
            # Explicitly turn off scientific notation.
            wigdf[wigdf[,1] == seqname,],
            file = out.file,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            dec = ".",
            append = TRUE
          );
        }
      }
      else {
        base::cat("  failed!\n");
      }
    }

    # dam format.
    if (DAM) {
      base::cat("creating dam file...\n");
      dam <- NULL;
      try(expr =
         dam <- GATCagg(
             MAnorm = MAnorm,  # Just created.
             marray = marray,
             targets = targets
         )
      );
      # Write wig data.frame to file.
      if (!is.null(dam)) {
        # Output file name.
        this.out.path <- switch(
            EXPR = sortby,
                "name" = file.path(out.path, protname),
                "filetype" = file.path(out.path, "dam"),
                "none" = out.path
        );
        if (!file.exists(this.out.path)) {
          dir.create(this.out.path);
        }
        out.file <- .escape(
             paste(core.name, "_", .dtag(), ".dam", sep=""));
        out.file <- file.path(this.out.path, out.file);
        # Add the 'meta' line to the vheader.
        addcomment(dam, "meta data", paste(meta[i,], collapse=" "));
        write.table(
          dam,
          file = out.file,
          row.names = FALSE,
          quote = FALSE,
          sep = "\t"
        );
      }
      else {
        base::cat("  failed!\n");
      }
    }


    ### various plots
    if (this.plotBias) {
      this.out.path <- switch(
          EXPR = sortby,
              "name" = file.path(out.path, protname),
              "filetype" = file.path(out.path, "plasmid-bias-plots"),
              "none" = out.path
      );
      if (!file.exists(this.out.path)) {
        dir.create(this.out.path);
      }
      base::cat("creating plasmid bias plot...\n");
      try(expr =
         plot.plasmid.bias(
             out.path = this.out.path,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             MAnorm = MAnorm,
             array1.name = array1.name,
             array2.name = array2.name,
             marray = marray,
             gene = gene,
             cex = cex,
             graph = graph
         )
      );
    }


    if (plotMA) {
      base::cat("creating MA plot...\n");
      this.out.path <- switch(
          EXPR = sortby,
              "name" = file.path(out.path, protname),
              "filetype" = file.path(out.path, "MA-plots"),
              "none" = out.path
      );
      if (!file.exists(this.out.path)) {
        dir.create(this.out.path);
      }
      try(expr =
         plot.MA(
             out.path = this.out.path,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             MAnorm = MAnorm,
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
      this.out.path <- switch(
          EXPR = sortby,
              "name" = file.path(out.path, protname),
              "filetype" = file.path(out.path, "scatter-plots"),
              "none" = out.path
      );
      if (!file.exists(this.out.path)) {
        dir.create(this.out.path);
      }
      try(
         scatterplot(
             out.path = this.out.path,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             MAnorm = MAnorm,
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
      this.out.path <- switch(
          EXPR = sortby,
              "name" = file.path(out.path, protname),
              "filetype" = file.path(out.path, "acf-plots"),
              "none" = out.path
      );
      if (!file.exists(this.out.path)) {
        dir.create(this.out.path);
      }
      try(expr =
         plot.acf.norm(
             out.path = this.out.path,
             marray = marray,
             core.name = core.name,
             name = protname,
             intensity = meta$intensity[i],
             MAnorm = MAnorm,
             array1.name = array1.name,
             array2.name = array2.name,
             maxlag = 200,
             cex = cex,
             graph = graph
        )
      );
    }

    if (plotSample) {
      base::cat("creating profiles...\n");
      this.out.path <- switch(
          EXPR = sortby,
              "name" = file.path(out.path, protname),
              "filetype" = file.path(out.path, "sample-plots"),
              "none" = out.path
      );
      if (!file.exists(this.out.path)) {
        dir.create(this.out.path);
      }
      try(expr =
         plot.profile(
             out.path = this.out.path,
             core.name = core.name,
             name = protname,
             MAnorm = MAnorm,
             intensity = meta$intensity[i],
             marray = marray,
             cex = cex,
             graph = graph
        )
      );
    }
  }
}
