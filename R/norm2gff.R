norm2gff <- function(name, core.name, MAnorm, marray) {

  # Remove non mappable probes.
  MAnorm <- MAnorm[
      complete.cases(MAnorm[,c("seqname", "start", "end")]),
  ];

  # Create the GFF and write to disk.
  gff <- data.frame(
      MAnorm$seqname,
      ".",
      paste(core.name, "_loess", sep=""),
      MAnorm$start,
      MAnorm$end,
      MAnorm$M.norm,
      ".",
      ".",
      MAnorm$probeID);

  # Version tracking (vtrackR).
  gff <- vtag(gff);
  addcomment(gff, "array platform", marray$name);
  addcomment(gff, "release", marray$release);

  return (gff);

}
