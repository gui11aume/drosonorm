GATCagg <- function (MAnorm, marray,
   GATCprobes = c("merge", "exclude", "keep")) {
# Average the scores of the probes mapping to the same GATC fragment.

# GATCprobes: specifies what to do with probes containing a GATC
#   merge: merge them in the leftmost fragment
#   exclude: discard the probes with a GATC
#   keep: Do not merge probes with a GATC in any fragment

  if (marray$release == 5) {
    lookup <- nimbleGenArrayGATCfragmentLookup.r5; # lazy loaded
    GATCmap <- DmelGATCfragmentsInfo.r5;           # lazy loaded
  }
  else if (marray$release == 4) {
    lookup <- nimbleGenArrayGATCfragmentLookup.r4; # lazy loaded
    GATCmap <- DmelGATCfragmentsInfo.r4;           # lazy loaded
  }
  else {
    stop (paste("release", marray$release, "not supported"));
  }

  GATCprobes = match.arg(arg = GATCprobes);

  # Get the indices of GATCprobes in the array specifications.
  GATCindex <- lookup$probeID %in%
      nimbleGenArrayGATCProbeIdentifiers;          # lazy loaded

  if (GATCprobes == "exclude") {
    # Set the associated GATC fragment to NA (discard the probe).
    lookup$GATCfragment[GATCindex] <- NA;
  }
  else if (GATCprobes == "keep") {
    # Set the associated GATC fragment to the GATC probe itself.
    lookup$GATCfragment[GATCindex] <- lookup$probeID[GATCindex];
  }
  else if (GATCprobes != "merge") {
    stop ('GATCprobes must be "exclude", "keep" or "merge"');
  }

  merged <- merge(MAnorm, lookup);
  avg <- tapply(
            X = merged$M.norm,
            INDEX = GATCfragment,
            FUN = mean,
            na.rm = TRUE
         );

  return(vtag(
      merge(
          GATCmap,
          data.frame(GATCfragment = names(agg), score = agg)
      )
  ));

}
