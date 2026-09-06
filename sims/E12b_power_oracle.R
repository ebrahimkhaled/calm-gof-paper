# E12b -- the two blocks E12 did not finish: (B) power by basis under the final degree rule, (D) oracle-Sigma_x check.
# Usage: Rscript E12b_power_oracle.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260910
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mkB <- function(name, n, p, lambda, ar = 0.7, norm = 2.31) { Sig <- ar^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)); list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, Sig = Sig, sd_eta = norm) }
mkA <- function(name) { b <- c(.5, -.4, .3, -.3, .2); list(name = name, n = 500, p = 5, lambda = 100, b0 = b, Sc = diag(5), Sig = diag(5), sd_eta = sqrt(sum(b^2))) }
JOBS <- list()
for (d in list(mkA("A_k01"), mkB("B_k10", 400, 40, 416), mkB("B_k25", 400, 100, 416)))
  for (dg in c(2, 3)) for (rl in c(1, 1.5)) JOBS[[length(JOBS) + 1]] <- list(d = d, deg = dg, rel = rl)
one <- function(job, r) {
  d <- job$d; set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc
  eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + d$sd_eta * job$rel * if (job$deg == 2) (Zs^2 - 1) / sqrt(2) else (Zs^3 - 3 * Zs) / sqrt(6)
  y <- rbinom(d$n, 1, cf_expit(eta))
  pc <- cf_pieces(X, y, d$lambda, G = 10); R5 <- cf_reference_R5(X, y, d$lambda, pc)$R5
  ad <- cf_adaptive_degree(X, y, pc); kappa <- d$p / d$n
  g <- function(nm) kappa * c(dec = 0.24, e1 = 0.11, e2 = 0.11, e3 = 1.15)[[nm]]
  pk <- function(k) { Zk <- cf_edge_basis(pc, k); cf_pvalue(cf_qf(pc$v, Zk), (1 + g(c("e1", "e2", "e3")[k])) * R5, Zk) }
  p2 <- pk(2); p3 <- pk(3)
  k_new <- if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L         # rho2k (final)
  k_old <- max(2L, ad$degree)                                       # rho2 (previous)
  data.frame(rep = r, cell = d$name, dep_deg = job$deg, rel = job$rel, rho_hat = ad$rho_hat, k_new = k_new, k_old = k_old,
             hl = cf_pvalue(pc$S_dec, (1 + g("dec")) * R5), e2 = p2, e3 = p3,
             edge_new = if (k_new == 3L) p3 else p2, edge_old = if (k_old == 3L) p3 else p2)
}
one_oracle <- function(d, r) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, cf_expit(drop(X %*% d$b0)))
  pc <- cf_pieces(X, y, d$lambda, G = 10)
  bh <- cf_bellec_adjustments(X, y, pc); bo <- cf_bellec_adjustments(X, y, pc, Sigma_x = d$Sig)
  data.frame(rep = r, cell = d$name, a_hat = sqrt(bh$a2), a_or = sqrt(bo$a2), s_hat = sqrt(bh$s2), s_or = sqrt(bo$s2),
             rho_hat = sqrt(bh$a2 / (bh$a2 + bh$s2)), rho_or = sqrt(bo$a2 / (bo$a2 + bo$s2)))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
pw <- do.call(rbind, lapply(JOBS, function(job) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(job, r), error = function(e) NULL)))
  cat(sprintf("power %s deg%d rel%.1f n=%s (%.1f min)\n", job$d$name, job$deg, job$rel, if (is.null(o)) "0" else nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
write.csv(pw, file.path(OUT, "E12b_power_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
pag <- aggregate(cbind(hl, e2, e3, edge_new, edge_old) ~ cell + dep_deg + rel, pw, rej)
pag$rho_hat <- round(aggregate(rho_hat ~ cell + dep_deg + rel, pw, mean)$rho_hat, 3)
pag$share3_new <- aggregate(k_new ~ cell + dep_deg + rel, pw, function(v) mean(v == 3))$k_new
pag$share3_old <- aggregate(k_old ~ cell + dep_deg + rel, pw, function(v) mean(v == 3))$k_old
pag$R <- R; write.csv(pag, file.path(OUT, "E12b_power_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== (B) power by basis, final rule vs previous | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(pag[order(pag$cell, pag$dep_deg, pag$rel), ], digits = 3, row.names = FALSE)
ORC <- list(mkB("B_k25", 400, 100, 416), mkB("B_k40", 400, 160, 416), mkB("H5_k25_ar03", 400, 100, 416, ar = 0.3))
orc <- do.call(rbind, lapply(ORC, function(d) do.call(rbind, parLapply(cl, seq_len(min(R, 200)), function(r) tryCatch(one_oracle(d, r), error = function(e) NULL)))))
stopCluster(cl); write.csv(orc, file.path(OUT, "E12b_oracle_perrep.csv"), row.names = FALSE)
oag <- aggregate(cbind(a_hat, a_or, s_hat, s_or, rho_hat, rho_or) ~ cell, orc, mean); write.csv(oag, file.path(OUT, "E12b_oracle_summary.csv"), row.names = FALSE)
cat("\n=== (D) sample-covariance plug-in vs oracle Sigma_x ===\n"); print(oag, digits = 4, row.names = FALSE)
