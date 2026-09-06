# E8 -- Package item D: the attenuation law as a DESIGN principle for the basis.
# Series C Prop 3: a degree-k departure survives grouping on the fitted index with non-centrality rho^{2k}.
# So at large kappa (low rho) a LOWER-degree basis should have more power.  rho is OBSERVABLE from Bellec's
# adjustments: rho_hat = a_hat / sqrt(a_hat^2 + sigma_hat^2)  (correlation of the fitted index with the true one).
# Bases on the corrected residual v (R2 reference, fitted grouping):  decile (all G directions), EDGE-1 (linear in
# group-mean pi_hat), EDGE-2 (deg <= 2), EDGE-3 (deg <= 3, the default).  Departures along the index: Hermite
# degree-2 and degree-3 at REL in {1.0, 1.5}; cells kappa .01 (A), .10 (p40), .25 (p100).  Size row REL = 0.
# Usage: Rscript E8_basis_degree_vs_kappa.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 500; W <- if (length(args) > 1) as.integer(args[2]) else 10
SEED <- 20260911
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0))) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416))
DEPS <- rbind(data.frame(deg = 0, rel = 0), expand.grid(deg = c(2, 3), rel = c(1.0, 1.5)))
He <- function(Z, k) switch(as.character(k), "2" = (Z^2 - 1) / sqrt(2), "3" = (Z^3 - 3 * Z) / sqrt(6))
one <- function(r, d, deg, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  if (deg > 0) eta <- eta + d$sd_eta * rel * He(Zs, deg)
  y <- rbinom(d$n, 1, cf_expit(eta))
  pc <- cf_pieces(X, y, d$lambda); rf <- cf_references(X, y, d$lambda, pc); Sg <- rf$R2
  b <- cf_bellec_adjustments(X, y, pc); rho_hat <- sqrt(b$a2 / (b$a2 + b$s2))
  pbar <- vapply(split(seq_len(d$n), pc$groups), function(I) mean(pc$pi[I]), 0.0)
  Zk <- lapply(1:3, function(k) as.matrix(stats::poly(pbar, k)))
  data.frame(rep = r, cell = d$name, deg = deg, rel = rel, rho_hat = rho_hat,
             p_dec = cf_pvalue(pc$S_dec, Sg),
             p_e1 = cf_pvalue(cf_qf(pc$v, Zk[[1]]), Sg, Zk[[1]]), p_e2 = cf_pvalue(cf_qf(pc$v, Zk[[2]]), Sg, Zk[[2]]), p_e3 = cf_pvalue(cf_qf(pc$v, Zk[[3]]), Sg, Zk[[3]]))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) do.call(rbind, parLapply(cl, seq_len(R), one, d = d, deg = DEPS$deg[k], rel = DEPS$rel[k]))))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E8_basis_degree_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(p_dec, p_e1, p_e2, p_e3) ~ cell + deg + rel, data = res, FUN = rej); agg$rho_hat <- aggregate(rho_hat ~ cell + deg + rel, data = res, FUN = median)$rho_hat
agg$R <- R; write.csv(agg, file.path(OUT, "E8_basis_degree_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E8 | basis degree vs kappa (R2 reference) | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("p_dec", "p_e1", "p_e2", "p_e3")) sh[[v]] <- sprintf("%.3f", sh[[v]]); sh$rho_hat <- sprintf("%.3f", sh$rho_hat)
print(sh[order(sh$cell, sh$deg, sh$rel), c("cell", "deg", "rel", "rho_hat", "p_dec", "p_e1", "p_e2", "p_e3")], row.names = FALSE)
cat("\nrho_hat = a/sqrt(a^2+sigma^2) (Bellec), median over reps: A ~.97, p40 ~.9, p100 ~.72 expected. rel = 0 rows are size. Prop 3 predicts EDGE-1/2 beat EDGE-3 as rho falls.\n")
