GATCagg <- function (MAnorm, marray,
   GATCprobes = c("merge", "exclude", "keep"), targets=TRUE) {
# Average the scores of the probes mapping to the same GATC fragment.

# GATCprobes: specifies what to do with probes containing a GATC
#   merge: merge them in the leftmost fragment
#   exclude: discard the probes with a GATC
#   keep: Do not merge probes with a GATC in any fragment

  if (marray$release == "dm3/R5") {
    lookup <- nimbleGenArrayGATCfragmentLookup.r5; # lazy loaded
    GATCmap <- DmelGATCfragmentsInfo.r5;           # lazy loaded
  }
  else if (marray$release == "dm2/R4") {
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
  # Rounding quickens read from disk tremendously.
  avg <- round(tapply(
            X = merged$M.norm,
            INDEX = merged$GATCfragment,
            FUN = mean,
            na.rm = TRUE
         ), 3);

  dam <- merge(GATCmap, data.frame(fragmentID=names(avg), score=avg));
  if(targets) {
    require(HMMt);
    bdg <- bridge(dam[,c("seqname", "start", "end", "score")]);
    fit <- BaumWelchT(
        x = bdg$x,
        series.length = bdg$series.length,
        m = 3
    );
    boundstate <- which.max(fit$mu);
    dam$bound <- 0 + (fit$ViterbiPath[bdg$nonvirtuals] == boundstate);
  }
  addcomment(dam, "array platform", marray$name);
  addcomment(dam, "release", marray$release);

  return (dam);

}
