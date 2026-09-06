# E3e -- Kill the SELECTION non-centrality by grouping on the one-step LEAVE-ONE-OUT fitted index.
# eta_loo_i = eta_hat_i - h_ii (y_i - pi_hat_i) / (w_i (1 - h_ii)),  h_ii = w_i x_i' M^{-1} x_i   (Pregibon 1981, penalized IRLS)
# The group of observation i no longer depends on its own y_i, so the self-selection that puts positive-noise
# observations into the top decile is removed to first order.  Everything else (r, U, mu_hat, Omega) as before,
# with the grouping held at the LOO deciles.  Deterministic, one fit, no resampling.
# Reports, per cell: deflation and size for R0 / R2 / R4 under LOO grouping, and ||mean v||^2 vs E3 (random grouping).
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 500; W <- if (length(args) > 1) as.integer(args[2]) else 6
SEED <- 20260905
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
design <- function(tag, n, p, lambda, rho, beta_norm) {
  if (tag == "A") { beta0 <- c(0.5, -0.4, 0.3, -0.3, 0.2); Sc <- diag(p) }
  else { beta0 <- rep_len(c(0.35, -0.30, 0.25, -0.20, 0.15), p); beta0 <- beta0 * beta_norm / sqrt(sum(beta0^2)); Sc <- chol(rho^abs(outer(1:p, 1:p, "-"))) }
  list(tag = tag, n = n, p = p, lambda = lambda, beta0 = beta0, Sc = Sc, name = sprintf("%s_n%d_p%d_lam%g", tag, n, p, lambda))
}
CELLS <- list(design("A", 500, 5, 100, 0, NA), design("A", 500, 5, 137, 0, NA), design("A", 500, 5, 200, 0, NA),
              design("B", 400, 100, 50, 0.7, 2.31), design("B", 400, 100, 416, 0.7, 2.31), design("B", 400, 100, 1000, 0.7, 2.31),
              design("B", 400, 40, 416, 0.7, 2.31), design("B", 400, 160, 416, 0.7, 2.31))
grp <- function(v, G = 10) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)
one_rep <- function(r, d) {
  set.seed(SEED + r)
  X <- matrix(rnorm(d$n * d$p), d$n, d$p) %*% d$Sc; eta0 <- drop(X %*% d$beta0); y <- rbinom(d$n, 1, cf_expit(eta0))
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, d$p)))
  beta <- cf_ridge(X1, y, d$lambda, Dpen); eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
  M <- crossprod(X1, w * X1) + d$lambda * Dpen
  h <- w * rowSums((X1 %*% solve(M)) * X1)                      # penalized hat-matrix diagonal
  eta_loo <- eta - h * (y - pi) / (w * pmax(1 - h, 1e-6))
  pc <- cf_pieces(X, y, d$lambda, beta = beta, groups = grp(eta_loo))
  ref <- cf_references(X, y, d$lambda, pc)
  out <- data.frame(rep = r, cell = d$name, S_dec = pc$S_dec, S_edge = pc$S_edge, mean_h = mean(h))
  for (k in c("R0", "R2", "R4")) { Sg <- ref[[k]]
    out[[paste0(k, "_df_dec")]] <- cf_implied_df(Sg); out[[paste0(k, "_df_edge")]] <- if (is.null(pc$Z)) NA else cf_implied_df(Sg, pc$Z)
    out[[paste0(k, "_p_dec")]] <- cf_pvalue(pc$S_dec, Sg); out[[paste0(k, "_p_edge")]] <- cf_pvalue(pc$S_edge, Sg, pc$Z) }
  vm <- matrix(pc$v, nrow = 1); colnames(vm) <- paste0("v", seq_along(pc$v)); cbind(out, vm)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one_rep, d = d)))); stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E3e_loo_perrep.csv"), row.names = FALSE)
rows <- list()
for (cn in unique(res$cell)) { s <- res[res$cell == cn, ]; vbar <- colMeans(s[, grep("^v[0-9]+$", names(s)), drop = FALSE])
  for (k in c("R0", "R2", "R4")) rows[[length(rows) + 1]] <- data.frame(cell = cn, ref = k,
    defl_dec = mean(s$S_dec) / mean(s[[paste0(k, "_df_dec")]]), size_dec = mean(s[[paste0(k, "_p_dec")]] < .05),
    defl_edge = mean(s$S_edge, na.rm = TRUE) / mean(s[[paste0(k, "_df_edge")]], na.rm = TRUE), size_edge = mean(s[[paste0(k, "_p_edge")]] < .05, na.rm = TRUE),
    mean_S = mean(s$S_dec), meanv2 = sum(vbar^2), mean_h = mean(s$mean_h)) }
summ <- do.call(rbind, rows); summ$R <- R; write.csv(summ, file.path(OUT, "E3e_loo_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E3e | LOO-index grouping | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- summ; for (v in c("defl_dec", "size_dec", "defl_edge", "size_edge", "meanv2", "mean_h")) sh[[v]] <- sprintf("%.3f", sh[[v]]); sh$mean_S <- sprintf("%.2f", sh$mean_S)
print(sh, row.names = FALSE)
cat("\nCompare meanv2 with E3 (random grouping): A ~.01, B50 .087, B416 .047*6.02=.28, B1000 .06*6.03=.36, p40 .024*7.25=.17, p160 .086*4.81=.41.  Right reference: defl 1.00, size .05.\n")
