.dtag <- function() {
# Date tag commodity function.
    format(Sys.Date(), "%y%m%d");
}

.checknames <- function(names) {
# Check that all needed column names are present
# in meta file, and 

  required <- c(
      "name",
      "exp1",
      "exp2",
      "ctl1",
      "ctl2",
      "exp1spike",
      "exp2spike",
      "ctl1spike",
      "ctl2spike",
      "intensity"
  )

  names <- tolower(names);
  missing <- which(! required %in% names);
  if (length(missing) > 0) {
    stop(paste(
          "missing column(s)",
          names[missing],
          "in meta information"
    ));
  }
  
  return (names);

}

.print.vcontrol <- function(pos=c("left", "right"), ...) {
    # Print version control in the top-left corner.
    pos <- match.arg(arg=pos);
    info <- c(
        paste("drosornom", packageDescription('drosonorm')$Revision),
        date()
    );
    if (pos == "left") {
      text <- paste(paste(" ", info, sep=""), collapse="\n");
      mtext(text, adj=0, line=-2, ...);
    }
    else if (pos == "right") {
      text <- paste(paste(info, " ", sep=""), collapse="\n");
      mtext(text, adj=1, line=-2, ...);
    }
  }

.escape <- function(string) {
  # Assert that the sring is valid UNIX file name.
  return (gsub("[][:space:]&%)(#`'\"!;:}{><,[]", "_", string));
}


.getgene <- function(gene, release) {
# Get gene mapping from gene ID.

  if (release == "dm3/R5") {
    genes <- genes.r5.17;                           # lazy loaded
  }
  else if (release == "dm2/R4") {
    genes <- genes.r4.3;                            # lazy loaded
  }
  else {
    stop ("non supported release");
  }

  # Identify the gene (lazy loaded).
  if (index <- match(gene, genes$geneID, nomatch = FALSE)) {
    return (list(
      seqname = genes$seqname[index],
      start = genes$start[index],
      end = genes$end[index]
    ));
  }
  else {
    return (NULL);
  }
}
