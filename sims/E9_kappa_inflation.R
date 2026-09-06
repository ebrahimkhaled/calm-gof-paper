# E9 -- "put kappa into the distribution": inflate the R5 reference covariance by (1 + g_hat), where g_hat is the
# predicted non-centrality fraction ||E v||^2 / tr(Omega) from OBSERVABLES (law fitted in sims/_delta_law.R on 8 cells):
#   K  : g = 0.24 * kappa                      (kappa = p/n)
#   H  : g = 0.28 * sqrt(hbar),  hbar = df_M/n  (mean self-influence)
#   HR : g = 0.50 * hbar + 0.07 * (1 - rho_hat^2)
# A scaled central weighted chi-square with the right mean has a HEAVIER upper tail than the true non-central law,
# so the adjustment can only err toward conservatism.  Size at 5% for R5 (unadjusted) and the three adjusted references.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 500; W <- if (length(args) > 1) as.integer(args[2]) else 6
SEED <- 20260905
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
design <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_n%d_p%d_lam%g", tag, n, p, lambda), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(design("A", 500, 5, 100), design("A", 500, 5, 137), design("A", 500, 5, 200), design("B", 400, 100, 50), design("B", 400, 100, 416), design("B", 400, 100, 1000), design("B", 400, 40, 416), design("B", 400, 160, 416))
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, plogis(drop(X %*% d$b0)))
  pc <- cf_pieces(X, y, d$lambda); r5 <- cf_reference_R5(X, y, d$lambda, pc); Sg <- r5$R5
  b <- cf_bellec_adjustments(X, y, pc); rho <- sqrt(b$a2 / (b$a2 + b$s2)); hbar <- unname(pc$obs["df_M"]) / d$n; kappa <- d$p / d$n
  g <- c(K = 0.24 * kappa, H = 0.28 * sqrt(hbar), HR = 0.50 * hbar + 0.07 * (1 - rho^2))
  out <- data.frame(rep = r, cell = d$name, kappa = kappa, hbar = hbar, rho = rho, p5_dec = cf_pvalue(pc$S_dec, Sg), p5_edge = cf_pvalue(pc$S_edge, Sg, pc$Z))
  for (k in names(g)) { out[[paste0("p", k, "_dec")]] <- cf_pvalue(pc$S_dec, (1 + g[k]) * Sg); out[[paste0("p", k, "_edge")]] <- cf_pvalue(pc$S_edge, (1 + g[k]) * Sg, pc$Z); out[[paste0("g", k)]] <- g[k] }
  out
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one, d = d)))); stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E9_kappa_inflation_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(p5_dec, p5_edge, pK_dec, pK_edge, pH_dec, pH_edge, pHR_dec, pHR_edge) ~ cell, data = res, FUN = rej)
agg$gK <- aggregate(gK ~ cell, res, mean)$gK; agg$gH <- aggregate(gH ~ cell, res, mean)$gH; agg$gHR <- aggregate(gHR ~ cell, res, mean)$gHR; agg$R <- R
write.csv(agg, file.path(OUT, "E9_kappa_inflation_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E9 | kappa-inflated R5 reference | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in names(sh)[grep("^p", names(sh))]) sh[[v]] <- sprintf("%.3f", sh[[v]]); for (v in c("gK", "gH", "gHR")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[, c("cell", "p5_dec", "pK_dec", "pH_dec", "pHR_dec", "p5_edge", "pK_edge", "pH_edge", "pHR_edge", "gK", "gH", "gHR")], row.names = FALSE)
cat("\np5 = R5 unadjusted; pK = (1+0.24 kappa) R5; pH = (1+0.28 sqrt(hbar)) R5; pHR = (1+0.50 hbar+0.07(1-rho^2)) R5.  Target .050; conservative (< .05) is acceptable, liberal is not.\n")
