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

