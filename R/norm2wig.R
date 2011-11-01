norm2wig <- function(MAnorm, marray) {
# Return a WIG data.frame from a given normalized data.frame.
# Alexey Pindyurin & Guillaume Filion

  # Calculate mean positions of the probes and store them
  # in the column "start".
  mid <- round((MAnorm$start + MAnorm$end)/2);

  # Remove unnecessary columns. For coordinates, leave only
  # "start" (in fact, "middle") positions.
  wigdf <- data.frame(
    seqname = MAnorm$seqname,
    position = mid,
    score = MAnorm$M.norm
  );

  # Remove NAs
  wigdf <- wigdf[complete.cases(wigdf)];

  # Version tracking (vtrackR).
  wigdf <- vtg(wigdf);
  addcomment(MAnorm, "array platform", marray$name);
  addcomment(MAnorm, "release", marray$release);

  return (wigdf);

}