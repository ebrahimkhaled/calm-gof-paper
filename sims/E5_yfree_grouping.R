# E5 -- Package item A (+B): remove the SELF-SELECTION non-centrality by grouping on an index that
# depends less (or not at all) on y, and add the deterministic second-order mean (delta method) as a
# Davies non-centrality.
# Groupings compared, all with the SAME statistic (r, U, mu_hat at the deployed lambda):
#   fitted   : deciles of eta_hat(lambda)                       (the test as used; self-influence h_ii = O(kappa))
#   heavyx5/x20/x100 : deciles of eta_hat(m * lambda)           (self-influence shrinks with m)
#   pc1      : deciles of the first principal component of X    (a function of X only -- zero selection)
#   oracleX  : deciles of eta(beta_hat(pi_0))  [simulation-only truth for the y-free ideal]
# References: R5 (de-noised variance) central, and R5 + delta2 (second-order mean, fixed-grouping delta method).
# Usage: Rscript E5_yfree_grouping.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 400; W <- if (length(args) > 1) as.integer(args[2]) else 10
SEED <- 20260905
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
design <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d_lam%g", tag, p, lambda), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(design("A", 500, 5, 100), design("B", 400, 100, 416), design("B", 400, 160, 416))
grp <- function(v, G = 10) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)
nc_pvalue <- function(S, Sigma, delta, Z = NULL) {   # non-central weighted chi-square tail
  e <- eigen(Sigma, symmetric = TRUE); keep <- e$values > 1e-9; Q <- e$vectors[, keep, drop = FALSE]; lam <- e$values[keep]
  Sh <- Q %*% (sqrt(lam) * t(Q)); Shinv <- Q %*% ((1 / sqrt(lam)) * t(Q))
  B <- if (is.null(Z)) Sh %*% Sh else { P <- Z %*% solve(crossprod(Z), t(Z)); Sh %*% P %*% Sh }
  eb <- eigen(B, symmetric = TRUE); kb <- eb$values > 1e-9
  d <- CompQuadForm::davies(S, lambda = eb$values[kb], delta = as.numeric(crossprod(eb$vectors[, kb, drop = FALSE], Shinv %*% delta))^2, lim = 1e5, acc = 1e-7)
  p <- d$Qq; if (d$ifault != 0 || is.na(p) || p < 0 || p > 1) p <- CompQuadForm::imhof(S, lambda = eb$values[kb], delta = as.numeric(crossprod(eb$vectors[, kb, drop = FALSE], Shinv %*% delta))^2)$Qq
  min(max(p, 0), 1) }
delta2 <- function(X1, pi_gen, lambda, Dpen, idx, h = 0.25) {   # second-order delta-method mean of v, grouping fixed
  bbar <- cf_ridge(X1, pi_gen, lambda, Dpen); pb <- plogis(drop(X1 %*% bbar)); Mi <- solve(crossprod(X1, pmax(pb * (1 - pb), 1e-8) * X1) + lambda * Dpen)
  f <- function(ys) cf_v_of_beta(ys, X1, bbar + drop(Mi %*% crossprod(X1, ys - pi_gen)), lambda, Dpen, idx)
  v0 <- f(pi_gen); acc <- 0; w <- pi_gen * (1 - pi_gen)
  for (i in seq_len(nrow(X1))) { e <- numeric(nrow(X1)); e[i] <- h; acc <- acc + w[i] * (f(pi_gen + e) - 2 * v0 + f(pi_gen - e)) / h^2 }
  v0 + 0.5 * acc }
one_rep <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta0 <- drop(X %*% d$b0); pi0 <- plogis(eta0); y <- rbinom(d$n, 1, pi0)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, d$p)))
  beta <- cf_ridge(X1, y, d$lambda, Dpen)
  pc1 <- drop(X %*% eigen(crossprod(scale(X, scale = FALSE)) / d$n, symmetric = TRUE)$vectors[, 1])
  groupings <- list(fitted = grp(drop(X1 %*% beta)),
                    heavyx5 = grp(drop(X1 %*% cf_ridge(X1, y, 5 * d$lambda, Dpen))),
                    heavyx20 = grp(drop(X1 %*% cf_ridge(X1, y, 20 * d$lambda, Dpen))),
                    heavyx100 = grp(drop(X1 %*% cf_ridge(X1, y, 100 * d$lambda, Dpen))),
                    pc1 = grp(pc1),
                    oracleX = grp(drop(X1 %*% cf_ridge(X1, pi0, d$lambda, Dpen))))
  rows <- list()
  for (gn in names(groupings)) {
    pc <- cf_pieces(X, y, d$lambda, beta = beta, groups = groupings[[gn]])
    r5 <- cf_reference_R5(X, y, d$lambda, pc); Sg <- r5$R5
    dn <- cf_denoised_w0(X, y, pc); b <- cf_bellec_adjustments(X, y, pc); k <- b$a2 / (b$a2 + b$s2); mz <- b$m * k; tau <- sqrt(b$a2 * b$s2 / (b$a2 + b$s2))
    gh <- statmod_gh(20); pi5 <- numeric(d$n); for (q in 1:20) pi5 <- pi5 + gh$w[q] * plogis(dn$alpha + dn$c * (mz + tau * gh$x[q])); pi5 <- pmin(pmax(pi5, 1e-6), 1 - 1e-6)
    dl <- tryCatch(delta2(X1, pi5, d$lambda, Dpen, split(seq_len(d$n), groupings[[gn]])), error = function(e) rep(NA, length(pc$v)))
    rows[[gn]] <- data.frame(rep = r, cell = d$name, grouping = gn, S_dec = pc$S_dec, S_edge = pc$S_edge,
      df_dec = cf_implied_df(Sg), df_edge = if (is.null(pc$Z)) NA else cf_implied_df(Sg, pc$Z),
      p_dec = cf_pvalue(pc$S_dec, Sg), p_edge = cf_pvalue(pc$S_edge, Sg, pc$Z),
      pnc_dec = if (anyNA(dl)) NA else nc_pvalue(pc$S_dec, Sg, dl), pnc_edge = if (anyNA(dl) || is.null(pc$Z)) NA else nc_pvalue(pc$S_edge, Sg, dl, pc$Z),
      delta2_norm2 = sum(dl^2), t(setNames(pc$v, paste0("v", seq_along(pc$v)))))
  }
  do.call(rbind, rows)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one_rep, d = d)))); stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E5_yfree_perrep.csv"), row.names = FALSE)
rows <- list()
for (cn in unique(res$cell)) for (gn in unique(res$grouping)) { s <- res[res$cell == cn & res$grouping == gn, ]; vbar <- colMeans(s[, grep("^v[0-9]+$", names(s)), drop = FALSE])
  rows[[length(rows) + 1]] <- data.frame(cell = cn, grouping = gn, defl_dec = mean(s$S_dec) / mean(s$df_dec), size_dec = mean(s$p_dec < .05), size_dec_nc = mean(s$pnc_dec < .05, na.rm = TRUE),
    defl_edge = mean(s$S_edge, na.rm = TRUE) / mean(s$df_edge, na.rm = TRUE), size_edge = mean(s$p_edge < .05, na.rm = TRUE), size_edge_nc = mean(s$pnc_edge < .05, na.rm = TRUE),
    meanv2 = sum(vbar^2), delta2 = mean(s$delta2_norm2, na.rm = TRUE), R = R) }
summ <- do.call(rbind, rows); write.csv(summ, file.path(OUT, "E5_yfree_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E5 | R5 reference | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- summ; for (v in c("defl_dec", "size_dec", "size_dec_nc", "defl_edge", "size_edge", "size_edge_nc", "meanv2", "delta2")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh, row.names = FALSE)
cat("\nsize_*_nc = R5 with the delta-method second-order mean as Davies non-centrality.  meanv2 = ||mean v||^2 (the real non-centrality);  delta2 = mean ||delta2||^2 (what the delta method predicts).  Target: size .05, meanv2 ~ .01-.08.\n")
