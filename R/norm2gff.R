norm2gff <- function(core.name, name, norm, marray) {

  # Output file name.
  out.file <- paste(core.name, "_", .dtag(), ".gff", sep="");

  # Map probes in "norm" data.frame
  m <- merge(x=norm, y=marray$mapping, by="PROBE_ID");
  m <- m[order(m$seqname, m$start),];

  # Make the GFF and write to disk.
  gff <- data.frame(
      m$seqname,
      ".",
      paste(core.name, "_loess", sep=""),
      m$start,
      m$end,
      m$M.norm,
      ".",
      ".",
      m$PROBE_ID);

  # Write and return gff.

  write.table(
      vtag(gff), # Version tracking.
      file = out.file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t",
  );

  return (vtag(gff));

}
