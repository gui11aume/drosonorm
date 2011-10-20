norm2gff <- function(name, core.name, MAnorm) {

  # Make the GFF and write to disk.
  gff <- data.frame(
      MAnorm$seqname,
      ".",
      paste(core.name, "_loess", sep=""),
      MAnorm$start,
      MAnorm$end,
      MAnorm$M.norm,
      ".",
      ".",
      MAnorm$PROBE_ID);

  return (vtag(gff));

}
