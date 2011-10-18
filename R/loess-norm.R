loess.norm <- function(core.name, out.path, exp1.file, ctl1.file,
    exp2.file='', ctl2.file='', marray) {
# Important: assume array names are of the form 'CMF_123456_...'


  # Output file name.
  out.file <- paste(core.name, "_norm_", .dtag(), ".txt", sep="");


  #################################################
  ##      normalize function for 2 channels      ##
  #################################################

  loess1array <- function(exp.file, ctl.file) {
  # Return a data.frame with probe IDs, single channel
  # values, A value, raw M and normalized M.

    # Read in the two .pair files.
    expt <- read.delim(exp.file, comment.char='#',
        colClasses=c(rep('character',4),rep('numeric',7)));
    ctrl <- read.delim(ctl.file, comment.char='#',
        colClasses=c(rep('character',4),rep('numeric',7)));

    # Compute M and A values.
    M <- log2(expt[,10] / ctrl[,10]);
    A <- .5 * log2(expt[,10] * ctrl[,10]);

    # Loess-normalize (approximate trace-hat.).
    M.norm <- loess(
        M ~ A,
        control = loess.control(trace.hat='approximate'),
        na.action = na.omit
    )$residuals;

    # NOTE: Rounding to 3 digits saves a lot of disk space.
    probes <- expt[,4];
    expt <- round(expt[,10], 3);
    ctrl <- round(ctrl[,10], 3);
    A <- round(A, 3);
    M <- round(M - median(M, na.rm=TRUE), 3);
    M.norm <- round(M.norm - median(M.norm, na.rm=TRUE), 3);

    return (data.frame(probes, expt, ctrl, A, M, M.norm));

  }


  # Normalize the first array.
  array1.name <- sub('CMF_([0-9]+)_.*','\\1', exp1.file);
  MA <- loess1array(exp1.file, ctl1.file);

  if (exp2.file == '') {
    # There is one array.
    names(MA) <- c("PROBE_ID", "exp1", "ctl1", "A", "M", "M.norm");
  }

  else {
    # There are two arrays.
    array2.name <- sub('CMF_([0-9]+)_.*','\\1', exp2.file);
    MA2 <- loess1array(exp2.file, ctl2.file);

    # Early return (with NA) if rownames differ.
    if (!identical(MA[,1], MA2[,1])) {
      warning(paste("rownames differ between files:", exp1.file,
          ctl1.file, exp2.file, ctl2.file, "(not processing)"));
      return (NA);
    }

    MA <- data.frame(MA[1:3], MA2[,2:3], MA[,4:6], MA2[,4:6]);
    names(MA) <- c("PROBE_ID", "exp1", "ctl1", "exp2", "ctl2",
        "A1", "M1", "M1.norm", "A2", "M2", "M2.norm");
    # Average M and A.
    MA$A <- round((MA$A1 + MA$A2)/2, 3);
    MA$M.norm <- round((MA$M1.norm + MA$M2.norm)/2, 3);
  }


  # Write the norm file (change directory and turn off the
  # warning triggered by 'commentline').
  MA <- vtag(MA);

  # Version tracking (the call contains the release).
  write.table(
      MA, 
      file = file.path(out.path, out.file),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
  );

  # Return the normalized file to the pipeline.
  return (MA);

}
