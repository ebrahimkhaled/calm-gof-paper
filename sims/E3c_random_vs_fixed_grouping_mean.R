# E3c -- is the non-centrality at kappa a RANDOM-GROUPING (selection) effect?
# Design B, lambda = 416, H0.  Per rep: v with (i) the test's random grouping on pi_hat(y);
# (ii) grouping FIXED at the deciles of the noise-free fit eta(beta_hat(pi_0)) -- a function of X only.
# If ||mean v||^2 collapses under (ii), the mean is a selection effect of grouping on y.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 300; W <- if (length(args) > 1) as.integer(args[2]) else 6
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
source(file.path(ROOT, "R", "closedform_engine.R"))
p <- 100; n <- 400; lambda <- 416; SEED <- 20260905
Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2))
grp <- function(v, G = 10) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)
one <- function(r) {
  set.seed(SEED + r); X <- matrix(rnorm(n * p), n) %*% Sc; eta0 <- drop(X %*% b0); pi0 <- plogis(eta0); y <- rbinom(n, 1, pi0)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, p)))
  bbar <- cf_ridge(X1, pi0, lambda, Dpen)                                  # noise-free fit: function of X only
  g_fixed <- grp(drop(X1 %*% bbar))
  pr <- cf_pieces(X, y, lambda); pf <- cf_pieces(X, y, lambda, groups = g_fixed)
  c(S_rand = pr$S_dec, S_fix = pf$S_dec, tr_rand = unname(pr$obs["tr_OmMLE"]), tr_fix = unname(pf$obs["tr_OmMLE"]), v_rand = pr$v, v_fix = pf$v)
}
cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, parLapply(cl, seq_len(R), one)); stopCluster(cl)
vr <- res[, grep("^v_rand", colnames(res))]; vf <- res[, grep("^v_fix", colnames(res))]
cat(sprintf("\n=== E3c | B lam416 | R = %d ===\nrandom grouping: mean S = %.3f, tr = %.3f, defl = %.3f, ||mean v||^2 = %.3f (noise floor ~ %.3f)\nfixed  grouping: mean S = %.3f, tr = %.3f, defl = %.3f, ||mean v||^2 = %.3f\n",
  R, mean(res[, "S_rand"]), mean(res[, "tr_rand"]), mean(res[, "S_rand"]) / mean(res[, "tr_rand"]), sum(colMeans(vr)^2), mean(res[, "tr_rand"]) / R,
  mean(res[, "S_fix"]), mean(res[, "tr_fix"]), mean(res[, "S_fix"]) / mean(res[, "tr_fix"]), sum(colMeans(vf)^2)))
cat("per-group mean v, random:", sprintf("%+.3f", colMeans(vr)), "\nper-group mean v, fixed :", sprintf("%+.3f", colMeans(vf)), "\n")
write.csv(as.data.frame(res), file.path(ROOT, "results", "E3c_grouping_mean_perrep.csv"), row.names = FALSE)
