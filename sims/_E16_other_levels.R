# The reference is validated at one point of the null distribution only: every size table in the
# paper and the supplement reports the rejection rate at 0.05. The applied section then reads
# p-values of 0.012 and 0.033 as magnitudes, which is a claim about the tail and not about a five
# per cent decision. The per-replication p-values behind Tables 2 and 3 are already in the shipped
# archive, so the check costs nothing but has never been run.
ROOT <- "C:/Users/ebrah/.gemini/Projects/PDFs/paper_seriesB_closedform/results"
rows <- list()
add <- function(f, cellcol, dcol, ecol, tag) {
  d <- read.csv(file.path(ROOT, f))
  for (cn in unique(d[[cellcol]])) {
    a <- d[d[[cellcol]] == cn, ]
    rows[[length(rows) + 1]] <<- data.frame(
      table = tag, cell = cn, R = nrow(a),
      D01 = mean(a[[dcol]] < .01, na.rm = TRUE), D05 = mean(a[[dcol]] < .05, na.rm = TRUE), D10 = mean(a[[dcol]] < .10, na.rm = TRUE),
      E01 = mean(a[[ecol]] < .01, na.rm = TRUE), E05 = mean(a[[ecol]] < .05, na.rm = TRUE), E10 = mean(a[[ecol]] < .10, na.rm = TRUE))
  }
}
add("E4c3_size_perrep.csv", "cell", "hl", "edge_k", "T2")     # the eight calibration cells
add("E10_heldout_perrep.csv", "cell", "hl", "edge", "T3")     # the seventeen held-out cells
d <- do.call(rbind, rows)
write.csv(d, file.path(ROOT, "E16_levels_summary.csv"), row.names = FALSE)
options(width = 200)
cat("\n=== rejection rate of a CORRECT model at three nominal levels ===\n")
cat("D = CALM-D (decile), E = CALM-E (adaptive EDGE). Nominal 0.01 / 0.05 / 0.10.\n\n")
print(d, digits = 3, row.names = FALSE)
cat("\n--- ranges over the twenty-five cells ---\n")
for (nm in c("D01", "D05", "D10", "E01", "E05", "E10"))
  cat(sprintf("  %-4s min %.3f  max %.3f  median %.3f\n", nm, min(d[[nm]]), max(d[[nm]]), median(d[[nm]])))
cat("\n--- the two tails worth naming ---\n")
cat(sprintf("  CALM-D worst at 0.01: %s = %.3f (%.1fx nominal); worst at 0.05: %s = %.3f (%.1fx)\n",
            d$cell[which.max(d$D01)], max(d$D01), max(d$D01) / .01,
            d$cell[which.max(d$D05)], max(d$D05), max(d$D05) / .05))
cat(sprintf("  CALM-E worst at 0.01: %.3f; cells with ZERO rejections at 0.01: %s\n",
            max(d$E01), paste(d$cell[d$E01 == 0], collapse = ", ")))
