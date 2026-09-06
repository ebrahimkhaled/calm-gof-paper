# E10 -- HELD-OUT validation of CALM (referee request): the inflation constants (0.24k, 1.15k, 0.11k) and the degree rule were
# calibrated on the eight Table-1 cells. Here CALM (final default) is run on cells NOT used for calibration:
#   H1 kappa = .15 (p = 60), H2 kappa = .33 (p = 132), H3 kappa = .40 lambda = 50, H4 kappa = .40 lambda = 1000,
#   H5 kappa = .25 with weak correlation (AR 0.3), H6 kappa = .25 with weak signal ||beta|| = 1.5, H7 kappa = .25 strong signal ||beta|| = 3.0,
#   H8 n = 800, p = 200 (kappa = .25, lambda = 832), H9 kappa = .25 G = 20 groups.
# Plus CALM under t3 and binary predictor laws in the three Table-2 cells (A, kappa .10, kappa .25) and a kappa .40 row (the referees:
# "Table 2 reports Sigma_5, not CALM"). Size only (null). R = 500, no bootstrap -> minutes.
# Also a SPARSE design for the GRP fairness point: kappa .25 with only 5 non-zero coefficients (GRP's home ground): CALM + GRP,
# null / quadratic rel 1 / offindex rel 1.5, R = 300 (GRP ~10 s per call).   Usage: Rscript E10_heldout_validation.R [R] [Rgrp] [workers]
suppressPackageStartupMessages({library(parallel); library(GRPtests)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; RG <- if (length(args) >= 2) as.integer(args[2]) else 300; W <- if (length(args) >= 3) as.integer(args[3]) else 8
SEED <- 20260910
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mkB <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, law = "gaussian", G = 10, sparse = FALSE) {
  Sig <- ar^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig)
  b0 <- if (sparse) c(rep_len(c(.35, -.3, .25, -.2, .15), 5), rep(0, p - 5)) else rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)); list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, law = law, G = G, sd_eta = norm) }
mkA <- function(name, law) list(name = name, n = 500, p = 5, lambda = 100, b0 = c(.5, -.4, .3, -.3, .2), Sc = diag(5), law = law, G = 10, sd_eta = sqrt(sum(c(.5, -.4, .3, -.3, .2)^2)))
CELLS <- list(mkB("H1_k15", 400, 60, 416), mkB("H2_k33", 400, 132, 416), mkB("H3_k40_lam50", 400, 160, 50), mkB("H4_k40_lam1000", 400, 160, 1000),
              mkB("H5_k25_ar03", 400, 100, 416, ar = 0.3), mkB("H6_k25_weak", 400, 100, 416, norm = 1.5), mkB("H7_k25_strong", 400, 100, 416, norm = 3.0),
              mkB("H8_n800_k25", 800, 200, 832), mkB("H9_k25_G20", 400, 100, 416, G = 20),
              mkA("A_t3", "t3"), mkA("A_binary", "binary"), mkB("k10_t3", 400, 40, 416, law = "t3"), mkB("k10_binary", 400, 40, 416, law = "binary"),
              mkB("k25_t3", 400, 100, 416, law = "t3"), mkB("k25_binary", 400, 100, 416, law = "binary"), mkB("k40_t3", 400, 160, 416, law = "t3"), mkB("k40_binary", 400, 160, 416, law = "binary"))
draw <- function(d, r) { set.seed(SEED + r)
  Z <- switch(d$law, gaussian = matrix(rnorm(d$n * d$p), d$n), t3 = matrix(rt(d$n * d$p, 3) / sqrt(3), d$n), binary = matrix((rbinom(d$n * d$p, 1, .3) - .3) / sqrt(.21), d$n))
  X <- Z %*% d$Sc; eta <- drop(X %*% d$b0); list(X = X, eta = eta) }
one <- function(r, d) { dd <- draw(d, r); y <- rbinom(d$n, 1, cf_expit(dd$eta)); cf <- calm.gof(dd$X, y, d$lambda, G = d$G)
  data.frame(rep = r, cell = d$name, hl = cf$SC.HL$p.value, edge = cf$SC.EDGE.adaptive$p.value, deg = cf$SC.EDGE.adaptive$degree, rho_hat = cf$SC.EDGE.adaptive$rho_hat, e3 = cf$SC.EDGE$p.value) }
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(GRPtests)})))
res <- do.call(rbind, lapply(CELLS, function(d) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d), error = function(e) NULL)))
  cat(sprintf("done %s n=%d (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
write.csv(res, file.path(OUT, "E10_heldout_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(hl, edge, e3) ~ cell, res, rej); agg$share_k3 <- aggregate(deg ~ cell, res, function(v) mean(v == 3))$deg; agg$rho_hat <- round(aggregate(rho_hat ~ cell, res, mean)$rho_hat, 2); agg$R <- R
agg <- agg[match(sapply(CELLS, `[[`, "name"), agg$cell), ]; write.csv(agg, file.path(OUT, "E10_heldout_summary.csv"), row.names = FALSE)
cat(sprintf("\n=== E10 | held-out size of the FINAL CALM | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins")))); print(agg, digits = 3, row.names = FALSE)
# ---- sparse design: GRP on home ground ----
S <- mkB("S_k25_sparse5", 400, 100, 416, sparse = TRUE)
DEPS <- rbind(data.frame(dep = "null", rel = 0), data.frame(dep = "quadratic", rel = 1), data.frame(dep = "offindex", rel = 1.5))
oneS <- function(r, dep, rel) { dd <- draw(S, r); Zs <- dd$eta / S$sd_eta
  eta <- dd$eta + switch(dep, null = 0, quadratic = S$sd_eta * rel * (Zs^2 - 1) / sqrt(2), offindex = (rel * S$sd_eta / sqrt(1 + .09)) * dd$X[, 1] * dd$X[, 2])
  y <- rbinom(S$n, 1, cf_expit(eta)); cf <- calm.gof(dd$X, y, S$lambda)
  grp <- tryCatch(as.numeric(GRPtest(dd$X, y, fam = "binomial")), error = function(e) NA_real_)
  data.frame(rep = r, cell = S$name, dep = dep, rel = rel, calm_hl = cf$SC.HL$p.value, calm_edge = cf$SC.EDGE.adaptive$p.value, deg = cf$SC.EDGE.adaptive$degree, rho_hat = cf$SC.EDGE.adaptive$rho_hat, grp = grp) }
clusterExport(cl, c("S", "DEPS", "oneS"))
rs <- do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) { o <- do.call(rbind, parLapply(cl, seq_len(RG), function(r) tryCatch(oneS(r, DEPS$dep[k], DEPS$rel[k]), error = function(e) NULL)))
  cat(sprintf("done sparse %s rel %.1f n=%d (%.1f min)\n", DEPS$dep[k], DEPS$rel[k], nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl); write.csv(rs, file.path(OUT, "E10_sparse_perrep.csv"), row.names = FALSE)
ag2 <- aggregate(cbind(calm_hl, calm_edge, grp) ~ cell + dep + rel, rs, rej); nul <- rs[rs$dep == "null", ]
ag2$grp_adj <- sapply(seq_len(nrow(ag2)), function(i) { cv <- quantile(nul$grp, .05, na.rm = TRUE); mean(rs$grp[rs$dep == ag2$dep[i] & rs$rel == ag2$rel[i]] <= cv, na.rm = TRUE) })
ag2$rho_hat <- round(aggregate(rho_hat ~ cell + dep + rel, rs, mean)$rho_hat, 2); ag2$R <- RG; write.csv(ag2, file.path(OUT, "E10_sparse_summary.csv"), row.names = FALSE)
cat("\n=== sparse design (5 non-zero of 100), kappa .25: CALM vs GRP ===\n"); print(ag2, digits = 3, row.names = FALSE)
