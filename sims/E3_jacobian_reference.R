# =====================================================================
# E3 -- Do observable, matrix-level references remove the deflation?
#
# Same fits and groupings as E2b (random grouping = the test as used).  Four
# reference covariances for v = r - mu_hat, factorial in two ingredients:
#   R0  A*  W(beta_hat)   first-order Jacobian, fitted variance     (= Omega_MLE; the Series C row (a))
#   R1  J   W(beta_hat)   EXACT Jacobian,       fitted variance
#   R2  A*  W(beta_tilde) first-order Jacobian, DEBIASED-generator variance
#   R3  J   W(beta_tilde) both
# For each: implied E[S] = tr(Sigma) (decile) / tr(Sigma^{1/2} P_Z Sigma^{1/2}) (EDGE),
# the deflation ratio mean(S)/implied, and the analytic size at 5% (Davies).
# A reference is RIGHT when its deflation is 1.00 and its size is 0.05, in
# every cell, including the kappa sweep.  Also recorded: ||mean(v)||^2 relative
# to tr(Omega_MLE), to check the problem is covariance, not mean.
#
# Usage: Rscript E3_jacobian_reference.R [R] [workers]
# =====================================================================
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R       <- if (length(args) >= 1) as.integer(args[1]) else 500
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 20
SEED <- 20260905                       # SAME stream as E2b: identical datasets
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))

design <- function(tag, n, p, lambda, rho, beta_norm) {
  if (tag == "A") { beta0 <- c(0.5, -0.4, 0.3, -0.3, 0.2); Sc <- diag(p) }
  else { beta0 <- rep_len(c(0.35, -0.30, 0.25, -0.20, 0.15), p); beta0 <- beta0 * beta_norm / sqrt(sum(beta0^2)); Sc <- chol(rho^abs(outer(1:p, 1:p, "-"))) }
  list(tag = tag, n = n, p = p, lambda = lambda, beta0 = beta0, Sc = Sc, name = sprintf("%s_n%d_p%d_lam%g", tag, n, p, lambda))
}
CELLS <- list(
  design("A", 500,   5, 100,  0,   NA), design("A", 500,   5, 137,  0,   NA), design("A", 500,   5, 200,  0,   NA),
  design("B", 400, 100,  50,  0.7, 2.31), design("B", 400, 100, 416, 0.7, 2.31), design("B", 400, 100, 1000, 0.7, 2.31),
  design("B", 400,  40, 416,  0.7, 2.31), design("B", 400, 160, 416, 0.7, 2.31))

one_rep <- function(r, d) {
  set.seed(SEED + r)
  X <- matrix(rnorm(d$n * d$p), d$n, d$p) %*% d$Sc
  eta0 <- drop(X %*% d$beta0); y <- rbinom(d$n, 1, cf_expit(eta0))
  pc  <- cf_pieces(X, y, d$lambda)
  ref <- cf_references(X, y, d$lambda, pc)
  out <- data.frame(rep = r, cell = d$name, S_dec = pc$S_dec, S_edge = pc$S_edge, w_ratio = ref$w_ratio, w_ratio_recal = ref$w_ratio_recal, recal_slope = ref$recal_slope,
                    v_norm2 = sum(pc$v^2), tr_OmMLE = pc$obs["tr_OmMLE"])
  for (k in c("R0", "R1", "R2", "R3", "R4")) {
    Sg <- ref[[k]]
    out[[paste0(k, "_df_dec")]]  <- cf_implied_df(Sg)
    out[[paste0(k, "_df_edge")]] <- if (is.null(pc$Z)) NA else cf_implied_df(Sg, pc$Z)
    out[[paste0(k, "_p_dec")]]   <- cf_pvalue(pc$S_dec, Sg)
    out[[paste0(k, "_p_edge")]]  <- cf_pvalue(pc$S_edge, Sg, pc$Z)
  }
  # the mean check needs v itself; store the G components
  vm <- matrix(pc$v, nrow = 1); colnames(vm) <- paste0("v", seq_along(pc$v))
  cbind(out, vm)
}

t0 <- Sys.time()
cl <- makeCluster(WORKERS); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one_rep, d = d))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E3_jacobian_perrep.csv"), row.names = FALSE)

# ---- summary ----------------------------------------------------------------
cells <- unique(res$cell); rows <- list()
for (cn in cells) {
  s <- res[res$cell == cn, ]
  vbar <- colMeans(s[, grep("^v\\d+$", names(s))])              # mean corrected residual per group
  for (k in c("R0", "R1", "R2", "R3", "R4")) {
    rows[[length(rows) + 1]] <- data.frame(cell = cn, ref = k,
      defl_dec  = mean(s$S_dec)  / mean(s[[paste0(k, "_df_dec")]]),
      defl_edge = mean(s$S_edge, na.rm = TRUE) / mean(s[[paste0(k, "_df_edge")]], na.rm = TRUE),
      size_dec  = mean(s[[paste0(k, "_p_dec")]]  < 0.05),
      size_edge = mean(s[[paste0(k, "_p_edge")]] < 0.05, na.rm = TRUE),
      implied_dec = mean(s[[paste0(k, "_df_dec")]]), mean_S_dec = mean(s$S_dec),
      w_ratio = mean(s$w_ratio), w_ratio_recal = mean(s$w_ratio_recal), recal_slope = mean(s$recal_slope), meanv_norm2_over_tr = sum(vbar^2) / mean(s$tr_OmMLE), R = R)
  }
}
summ <- do.call(rbind, rows)
write.csv(summ, file.path(OUT, "E3_jacobian_summary.csv"), row.names = FALSE)

options(width = 200)
cat(sprintf("\n=== E3 | R = %d per cell | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
show <- summ
for (v in c("defl_dec", "defl_edge", "size_dec", "size_edge", "w_ratio", "w_ratio_recal", "recal_slope", "meanv_norm2_over_tr")) show[[v]] <- sprintf("%.3f", show[[v]])
for (v in c("implied_dec", "mean_S_dec")) show[[v]] <- sprintf("%.2f", show[[v]])
print(show[, c("cell", "ref", "defl_dec", "size_dec", "defl_edge", "size_edge", "implied_dec", "mean_S_dec", "w_ratio", "w_ratio_recal", "recal_slope", "meanv_norm2_over_tr")], row.names = FALSE)
cat("\nR0 = A* W(b_hat) [= Omega_MLE], R1 = J W(b_hat), R2 = A* W(b_tilde), R3 = J W(b_tilde), R4 = A* W(recalibrated scale a + b*eta_hat).",
    "\nA reference is right when defl = 1.00 and size = 0.05.  w_ratio = sum w(b_tilde)/sum w(b_hat).",
    "\nmeanv_norm2_over_tr = ||mean(v)||^2 / tr(Omega_MLE): ~0 means the corrected residual is centred (covariance problem only).\n")
