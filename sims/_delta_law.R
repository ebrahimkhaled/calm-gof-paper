# Does the selection non-centrality ||E v||^2 / tr(Omega) follow a law in OBSERVABLES across cells?
e3  <- read.csv("results/E3_jacobian_perrep_R0toR3_run1.csv")     # v columns, random grouping (the test)
e3b <- read.csv("results/E3b_denoised_perrep.csv")                  # a2/s2 not stored -- rho from est_signal_sd? no; recompute below
e2b <- read.csv("results/E2b_deflation_summary.csv"); e2b <- e2b[e2b$grouping == "random", ]
rows <- list()
for (cn in unique(e3$cell)) {
  s <- e3[e3$cell == cn, ]; vbar <- colMeans(s[, grep("^v[0-9]+$", names(s))]); tr <- mean(s$tr_OmMLE)
  g <- e2b[e2b$cell == cn, ]
  parts <- strsplit(cn, "_")[[1]]; n <- as.integer(sub("n", "", parts[2])); p <- as.integer(sub("p", "", parts[3])); lam <- as.numeric(sub("lam", "", parts[4]))
  rows[[length(rows) + 1]] <- data.frame(cell = cn, n = n, p = p, kappa = p / n, lambda = lam, meanv2_tr = sum(vbar^2) / tr,
    hbar = g$df_M / n, df_M = g$df_M, shrink = g$norm_bhat / g$norm_btilde, tr = tr)
}
d <- do.call(rbind, rows)
# rho_hat per cell: recompute from 60 draws per cell (cheap) via the engine
source("R/closedform_engine.R")
rho_cell <- function(n, p, lam) { Sc <- if (p == 5) diag(p) else chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- if (p == 5) c(.5, -.4, .3, -.3, .2) else { b <- rep_len(c(.35, -.3, .25, -.2, .15), p); b * 2.31 / sqrt(sum(b^2)) }
  median(replicate(60, { X <- matrix(rnorm(n * p), n) %*% Sc; y <- rbinom(n, 1, plogis(drop(X %*% b0))); pc <- cf_pieces(X, y, lam); b <- cf_bellec_adjustments(X, y, pc); sqrt(b$a2 / (b$a2 + b$s2)) })) }
set.seed(7); d$rho <- mapply(rho_cell, d$n, d$p, d$lambda)
d$one_minus_rho2 <- 1 - d$rho^2
options(width = 200); print(d[order(d$kappa, d$lambda), c("cell", "kappa", "lambda", "meanv2_tr", "hbar", "rho", "one_minus_rho2", "shrink")], row.names = FALSE, digits = 3)
cat("\ncandidate laws (fit through origin), R^2:\n")
for (f in c("hbar", "one_minus_rho2", "I(hbar*one_minus_rho2)", "I(one_minus_rho2^2)", "I(hbar^0.5)", "I(kappa)")) {
  m <- lm(as.formula(paste("meanv2_tr ~ 0 +", f)), data = d); cat(sprintf("  %-26s coef %.3f   R2 %.3f   max|resid| %.3f\n", f, coef(m), 1 - sum(resid(m)^2) / sum(d$meanv2_tr^2), max(abs(resid(m))))) }
m2 <- lm(meanv2_tr ~ 0 + hbar + one_minus_rho2, data = d); cat(sprintf("  hbar + (1-rho^2)          coefs %.3f %.3f  R2 %.3f\n", coef(m2)[1], coef(m2)[2], 1 - sum(resid(m2)^2) / sum(d$meanv2_tr^2)))
write.csv(d, "results/delta_law_cells.csv", row.names = FALSE)
