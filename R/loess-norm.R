loess.norm <- function(exp1.file, ctl1.file, exp2.file='',
    ctl2.file='', marray) {
# Important: assume array names are of the form 'CMF_123456_...'


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
    names(MA) <- c("probeID", "exp1", "ctl1", "A", "M", "M.norm");
  }

  else {
    # There are two arrays.
    array2.name <- sub('CMF_([0-9]+)_.*','\\1', exp2.file);
    MA2 <- loess1array(exp2.file, ctl2.file);

    # Sort on probe names (will sort on position later).
    MA <- MA[order(MA$probes),];
    MA2 <- MA2[order(MA2$probes),];

    # Early return (with NA) if rownames differ.
    if (!identical(MA[,1], MA2[,1])) {
      warning(paste("\nrownames differ between files:", exp1.file,
          ctl1.file, exp2.file, ctl2.file, "(not processing)\n"));
      return (NA);
    }

    MA <- data.frame(MA[1:3], MA2[,2:3], MA[,4:6], MA2[,4:6]);
    names(MA) <- c("probeID", "exp1", "ctl1", "exp2", "ctl2",
        "A1", "M1", "M1.norm", "A2", "M2", "M2.norm");
    # Average M and A.
    MA$A <- round((MA$A1 + MA$A2)/2, 3);
    MA$M.norm <- round((MA$M1.norm + MA$M2.norm)/2, 3);
  }

  # Save spikes and random probes in extraprobes for later.
  row.names(MA) <- MA$probeID;
  extraprobes <- MA[c(marray$spikes, marray$random),];
  extraprobes$seqname <- NA;
  extraprobes$start <- NA;
  extraprobes$end <- NA;

  # Map the probes in the given release and sort.
  # NB: this removes the double density, spikes and random probes...
  MAnorm <- merge(x=MA, y=marray$mapping, by="probeID");
  MAnorm <- MAnorm[order(MAnorm$seqname, MAnorm$start),];
  # ... so we add random probes and spikes back.
  MAnorm <- rbind(MAnorm, extraprobes);

  # Version tracking (vtrackR).
  MAnorm <- vtag(MAnorm);
  addcomment(MAnorm, "array platform", marray$name);
  addcomment(MAnorm, "release", marray$release);

  # Return the normalized file to the pipeline.
  return (MAnorm);

}
