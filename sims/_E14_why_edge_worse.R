# Post-processing of E14a (R = 500, already written): WHY does the smooth statistic fail harder
# than the decile one at low prevalence, when both see the same selection mean m = E(v)?
#
# The arithmetic that should decide it.  Writing tr for the reference trace and P for the projection
# onto the degree-k basis, the two non-centrality shares are ||m||^2 / tr(Sigma_0) and
# ||P m||^2 / tr(P Sigma_0).  Their ratio factors into
#     (fraction of ||m||^2 that lies IN the basis)  x  (tr Sigma_0 / tr(P Sigma_0)),
# and the second factor is about G/k because the projection keeps k of G directions.  So if the
# selection mean is a SMOOTH function of the group index -- which is what the self-influence argument
# predicts, and what E9/E9b found -- the first factor is near one, the smooth statistic inherits the
# whole displacement on a fifth of the trace, and it is worse by construction, at every prevalence.
# The balanced case hides this only because ||m||^2 is small there in absolute terms.
ROOT <- "C:/Users/ebrah/.gemini/Projects/PDFs/paper_seriesB_closedform"
A <- read.csv(file.path(ROOT, "results", "E14a_eta_basis_perrep.csv"))
G <- 10
out <- do.call(rbind, lapply(split(A, A$cell), function(D) {
  V <- as.matrix(D[, paste0("v", 1:G)]); m <- colMeans(V); k <- D$deg[1]
  fr <- mean(vapply(seq_len(nrow(D)), function(i) {
    Z <- as.matrix(stats::poly(as.numeric(D[i, paste0("pb", 1:G)]), k))
    drop(crossprod(crossprod(Z, m), solve(crossprod(Z), crossprod(Z, m)))) }, 0.0)) / sum(m^2)
  data.frame(cell = D$cell[1], prev = round(mean(D$prev), 3), k = k,
             norm_m2 = sum(m^2),                       # size of the selection displacement
             frac_in_basis = fr,                       # how much of it the smooth basis keeps
             trace_ratio = mean(D$tr_p) / (mean(D$tr_p) * 0 + mean(sapply(seq_len(nrow(D)), function(i) 1))) )
}))
# trace inflation factor: tr(Sigma_0) vs tr(P Sigma_0).  tr(Sigma_0) is not stored per rep, so take it
# from the decile ratio identity: ratio_d = mean(S_d)/tr(Sigma_0) is in the E14b/E14c summaries; here
# report the directly useful quantity instead -- the predicted and observed ratio of the two shares.
S <- read.csv(file.path(ROOT, "results", "E14a_eta_basis_summary.csv"))
out$ncsh_edge <- S$ncsh_p[match(out$cell, S$cell)]
out$pred_x <- out$frac_in_basis * (G / out$k)      # predicted (edge share) / (decile share)
cat("\n=== why the smooth statistic is hit harder: the selection mean is smooth, and the basis keeps it on a fifth of the trace ===\n\n")
print(out[order(-out$prev), c("cell", "prev", "k", "norm_m2", "frac_in_basis", "ncsh_edge", "pred_x")],
      digits = 3, row.names = FALSE)
cat("\nfrac_in_basis near 1 means the displacement is essentially a low-order polynomial in the group index.\n")
cat("pred_x = frac_in_basis * G/k is the factor by which the smooth statistic's non-centrality share should exceed the decile one.\n")
cat("\n--- the shape of the selection mean m = E(v), by group, low group first ---\n")
for (nm in c("k25_prev50", "k25_prev15", "k25_prev08")) {
  D <- A[A$cell == nm, ]; if (!nrow(D)) next
  cat(sprintf("%-12s %s\n", nm, paste(sprintf("%6.2f", colMeans(as.matrix(D[, paste0("v", 1:G)]))), collapse = " ")))
}
