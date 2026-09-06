# E10b -- diagnose the held-out failures of E10: H5 (kappa .25, AR 0.3: HL size .69!), H8 (n = 800, p = 200: .12), H3 (kappa .40, lambda 50: .13),
# against the calibration cell B kappa .25 (AR 0.7).  Per rep: S_dec, traces of R0/R2/R5, ||mean v||, Bellec (a, sigma, c), sd(eta_hat/eta_tilde/eta_0),
# p-values under R0 / R2 / R5 / CALM.  Plus prepivot (shrink.gof, B = 199) on the first 100 reps of H5 and H8 to see whether the bootstrap fails too.
suppressPackageStartupMessages({library(parallel); library(ebrahim.gof)})
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args) >= 1) as.integer(args[1]) else 200; RP <- if (length(args) >= 2) as.integer(args[2]) else 100; W <- 8
SEED <- 20260910
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mkB <- function(name, n, p, lambda, ar = 0.7, norm = 2.31) { Sig <- ar^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)); list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(mkB("B_k25_ar07", 400, 100, 416), mkB("H5_k25_ar03", 400, 100, 416, ar = 0.3), mkB("H8_n800_k25", 800, 200, 832), mkB("H3_k40_lam50", 400, 160, 50))
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta0 <- drop(X %*% d$b0); y <- rbinom(d$n, 1, cf_expit(eta0))
  pc <- cf_pieces(X, y, d$lambda, G = 10); R0 <- pc$Om_MLE; R2 <- cf_references(X, y, d$lambda, pc)$R2; r5 <- cf_reference_R5(X, y, d$lambda, pc); R5 <- r5$R5
  b <- cf_bellec_adjustments(X, y, pc); X1 <- cbind(1, X); D <- diag(c(0, rep(1, d$p))); bh <- pc$beta; eh <- drop(X1 %*% bh)
  wh <- cf_expit(eh) * (1 - cf_expit(eh)); Fm <- crossprod(X1, wh * X1); bt <- bh + solve(Fm, d$lambda * D %*% bh); et <- drop(X1 %*% bt)
  kappa <- d$p / d$n
  data.frame(rep = r, cell = d$name, S_dec = pc$S_dec, tr0 = sum(diag(R0)), tr2 = sum(diag(R2)), tr5 = sum(diag(R5)), v2 = sum(pc$v^2), meanv = paste(round(pc$v, 2), collapse = " "),
             a_hat = sqrt(b$a2), s_hat = sqrt(b$s2), rho_hat = sqrt(b$a2 / (b$a2 + b$s2)), c_hat = r5$c, sd_eta0 = sd(eta0), sd_eh = sd(eh), sd_et = sd(et), est_signal_sd = r5$est_signal_sd,
             norm_bh = sqrt(sum(bh[-1]^2)), norm_bt = sqrt(sum(bt[-1]^2)), df_M = unname(pc$obs["df_M"]),
             p0 = cf_pvalue(pc$S_dec, R0), p2 = cf_pvalue(pc$S_dec, R2), p5 = cf_pvalue(pc$S_dec, R5), pcalm = cf_pvalue(pc$S_dec, (1 + .24 * kappa) * R5))
}
onep <- function(r, d) { set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, cf_expit(drop(X %*% d$b0)))
  pp <- shrink.gof(X, y, lambda = d$lambda, G = 10, basis = "decile", B = 199, seed = SEED + 1e6 + r); data.frame(rep = r, cell = d$name, pp_hl = pp$SC.HL$p.value) }
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(ebrahim.gof)})))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d), error = function(e) NULL)))))
write.csv(res, file.path(OUT, "E10b_diag_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
options(width = 250)
agg <- aggregate(cbind(S_dec, tr0, tr2, tr5, v2, a_hat, s_hat, rho_hat, c_hat, sd_eta0, sd_eh, sd_et, est_signal_sd, norm_bh, norm_bt, df_M) ~ cell, res, mean)
sz <- aggregate(cbind(p0, p2, p5, pcalm) ~ cell, res, rej); agg <- merge(agg, sz, by = "cell")
cat(sprintf("\n=== E10b | diagnostics | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins")))); print(agg, digits = 3, row.names = FALSE)
cat("\nmean corrected residual by decile (first cell rows):\n"); for (c in unique(res$cell)) { vv <- do.call(rbind, lapply(strsplit(res$meanv[res$cell == c], " "), as.numeric)); cat(sprintf("%-14s", c), round(colMeans(vv), 3), "\n") }
pp <- do.call(rbind, lapply(CELLS[2:3], function(d) do.call(rbind, parLapply(cl, seq_len(RP), function(r) tryCatch(onep(r, d), error = function(e) NULL)))))
stopCluster(cl); write.csv(pp, file.path(OUT, "E10b_prepivot_perrep.csv"), row.names = FALSE)
cat("\nprepivot SC.HL size (B = 199):\n"); print(aggregate(pp_hl ~ cell, pp, rej), row.names = FALSE)
write.csv(agg, file.path(OUT, "E10b_diag_summary.csv"), row.names = FALSE)
