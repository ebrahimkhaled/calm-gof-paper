# =====================================================================
# E3b -- the de-noised variance (R5) and the second-order NON-CENTRALITY, together.
#
# E3 showed: R2 = A* W(beta_tilde) A*' is exact at fixed p and kappa = 0.10 and too
# TIGHT at kappa >= 0.25 (size .07-.13), and that at kappa > 0 the corrected residual
# v = r - mu_hat has a small non-zero MEAN (||mean v||^2 / tr = .05-.09).  So the
# reference at kappa needs (i) the right variance and (ii) a non-centrality.
#
# Variances tested:   R2  W at pi(beta_tilde)          R4  W at the recalibrated scale
#                     R5  W0 = E[pi0(1-pi0) | fit], de-noised via Bellec Thm 4.3 (engine)
# Non-centrality:     delta = v evaluated on the NOISE-FREE response y* = pi*(generator),
#                     i.e. the whole pipeline (IRLS, mu_hat) run on a real-valued y* in (0,1):
#                     this is the deterministic second-order displacement remainder that the
#                     bootstrap reproduces and the Taylor expansion overshot.  Two generators:
#                     y* = pi(beta_tilde)  and  y* = the R5 de-noised pi0-hat.
# References:         R2, R4, R5 (central) and R2nc, R5nc (non-central, Davies with delta).
# Same seeds and cells as E2b/E3.   Usage: Rscript E3b_denoised_noncentral.R [R] [workers]
# =====================================================================
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R       <- if (length(args) >= 1) as.integer(args[1]) else 500
WORKERS <- if (length(args) >= 2) as.integer(args[2]) else 20
SEED <- 20260905
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

# non-central weighted chi-square tail: S = u'Au, u ~ N(delta, Sigma), A = I or P_Z
cf_pvalue_nc <- function(S, Sigma, delta, Z = NULL) {
  if (is.na(S)) return(NA_real_)
  e <- eigen(Sigma, symmetric = TRUE); keep <- e$values > 1e-9
  Q <- e$vectors[, keep, drop = FALSE]; lam <- e$values[keep]
  Sh <- Q %*% (sqrt(lam) * t(Q)); Shinv <- Q %*% ((1 / sqrt(lam)) * t(Q))
  if (is.null(Z)) { B <- Sh %*% Sh } else { P <- Z %*% solve(crossprod(Z), t(Z)); B <- Sh %*% P %*% Sh }
  eb <- eigen(B, symmetric = TRUE); kb <- eb$values > 1e-9
  w  <- eb$values[kb]; nc <- as.numeric(crossprod(eb$vectors[, kb, drop = FALSE], Shinv %*% delta))^2
  d <- CompQuadForm::davies(S, lambda = w, delta = nc, lim = 1e5, acc = 1e-7); p <- d$Qq
  if (d$ifault != 0 || is.na(p) || p < 0 || p > 1) p <- CompQuadForm::imhof(S, lambda = w, delta = nc)$Qq
  min(max(p, 0), 1)
}

# v evaluated on a real-valued response y* (the noise-free pipeline), grouping held at the observed one
v_deterministic <- function(X, ystar, lambda, groups) {
  pcs <- cf_pieces(X, ystar, lambda, groups = groups)      # IRLS on real y* in (0,1) is well defined
  pcs$v
}

one_rep <- function(r, d) {
  set.seed(SEED + r)
  X <- matrix(rnorm(d$n * d$p), d$n, d$p) %*% d$Sc
  eta0 <- drop(X %*% d$beta0); y <- rbinom(d$n, 1, cf_expit(eta0))
  pc  <- cf_pieces(X, y, d$lambda)
  ref <- cf_references(X, y, d$lambda, pc)
  r5  <- cf_reference_R5(X, y, d$lambda, pc)
  X1  <- cbind(1, X)
  # generators for the non-centrality
  pi_t <- pmin(pmax(cf_expit(drop(X1 %*% pc$beta_tilde)), 1e-6), 1 - 1e-6)
  dn   <- cf_denoised_w0(X, y, pc)
  b    <- cf_bellec_adjustments(X, y, pc); k <- b$a2 / (b$a2 + b$s2); mu_z <- b$m * k; tau <- sqrt(b$a2 * b$s2 / (b$a2 + b$s2))
  gh <- statmod_gh(20); pi0_hat <- numeric(d$n)
  for (q in 1:20) pi0_hat <- pi0_hat + gh$w[q] * cf_expit(dn$alpha + dn$c * (mu_z + tau * gh$x[q]))
  pi0_hat <- pmin(pmax(pi0_hat, 1e-6), 1 - 1e-6)
  delta_t <- tryCatch(v_deterministic(X, pi_t,    d$lambda, pc$groups), error = function(e) rep(NA, length(pc$v)))
  delta_5 <- tryCatch(v_deterministic(X, pi0_hat, d$lambda, pc$groups), error = function(e) rep(NA, length(pc$v)))
  out <- data.frame(rep = r, cell = d$name, S_dec = pc$S_dec, S_edge = pc$S_edge, tr_OmMLE = pc$obs["tr_OmMLE"],
                    true_sd_eta = sd(eta0), est_signal_sd = r5$est_signal_sd, w_ratio_R2 = ref$w_ratio, w_ratio_R4 = ref$w_ratio_recal, w_ratio_R5 = r5$w_ratio_R5,
                    delta_t_norm2 = sum(delta_t^2), delta_5_norm2 = sum(delta_5^2))
  refs <- list(R2 = ref$R2, R4 = ref$R4, R5 = r5$R5)
  for (k2 in names(refs)) {
    Sg <- refs[[k2]]
    out[[paste0(k2, "_df_dec")]] <- cf_implied_df(Sg); out[[paste0(k2, "_df_edge")]] <- if (is.null(pc$Z)) NA else cf_implied_df(Sg, pc$Z)
    out[[paste0(k2, "_p_dec")]] <- cf_pvalue(pc$S_dec, Sg); out[[paste0(k2, "_p_edge")]] <- cf_pvalue(pc$S_edge, Sg, pc$Z)
  }
  # non-central variants: variance R2 with delta_t ; variance R5 with delta_5
  out$R2nc_p_dec  <- if (anyNA(delta_t)) NA else cf_pvalue_nc(pc$S_dec, ref$R2, delta_t)
  out$R2nc_p_edge <- if (anyNA(delta_t) || is.null(pc$Z)) NA else cf_pvalue_nc(pc$S_edge, ref$R2, delta_t, pc$Z)
  out$R5nc_p_dec  <- if (anyNA(delta_5)) NA else cf_pvalue_nc(pc$S_dec, r5$R5, delta_5)
  out$R5nc_p_edge <- if (anyNA(delta_5) || is.null(pc$Z)) NA else cf_pvalue_nc(pc$S_edge, r5$R5, delta_5, pc$Z)
  vm <- matrix(pc$v, nrow = 1); colnames(vm) <- paste0("v", seq_along(pc$v)); cbind(out, vm)
}

t0 <- Sys.time()
cl <- makeCluster(WORKERS); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one_rep, d = d))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E3b_denoised_perrep.csv"), row.names = FALSE)

rows <- list()
for (cn in unique(res$cell)) {
  s <- res[res$cell == cn, ]; vbar <- colMeans(s[, grep("^v[0-9]+$", names(s)), drop = FALSE])
  base <- data.frame(cell = cn, mean_S = mean(s$S_dec), tr_OmMLE = mean(s$tr_OmMLE), meanv2 = sum(vbar^2),
                     delta_t2 = mean(s$delta_t_norm2, na.rm = TRUE), delta_52 = mean(s$delta_5_norm2, na.rm = TRUE),
                     true_sd = mean(s$true_sd_eta), est_sd = mean(s$est_signal_sd),
                     wR2 = mean(s$w_ratio_R2), wR4 = mean(s$w_ratio_R4), wR5 = mean(s$w_ratio_R5))
  for (k2 in c("R2", "R4", "R5")) rows[[length(rows) + 1]] <- cbind(base, ref = k2,
      defl_dec = mean(s$S_dec) / mean(s[[paste0(k2, "_df_dec")]]), size_dec = mean(s[[paste0(k2, "_p_dec")]] < 0.05),
      defl_edge = mean(s$S_edge, na.rm = TRUE) / mean(s[[paste0(k2, "_df_edge")]], na.rm = TRUE), size_edge = mean(s[[paste0(k2, "_p_edge")]] < 0.05, na.rm = TRUE))
  for (k2 in c("R2nc", "R5nc")) rows[[length(rows) + 1]] <- cbind(base, ref = k2, defl_dec = NA, size_dec = mean(s[[paste0(k2, "_p_dec")]] < 0.05, na.rm = TRUE),
      defl_edge = NA, size_edge = mean(s[[paste0(k2, "_p_edge")]] < 0.05, na.rm = TRUE))
}
summ <- do.call(rbind, rows); summ$R <- R
write.csv(summ, file.path(OUT, "E3b_denoised_summary.csv"), row.names = FALSE)

options(width = 220)
cat(sprintf("\n=== E3b | R = %d per cell | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
show <- summ[, c("cell", "ref", "defl_dec", "size_dec", "defl_edge", "size_edge", "mean_S", "meanv2", "delta_t2", "delta_52", "true_sd", "est_sd", "wR2", "wR4", "wR5")]
for (v in c("defl_dec", "size_dec", "defl_edge", "size_edge", "wR2", "wR4", "wR5")) show[[v]] <- ifelse(is.na(show[[v]]), "", sprintf("%.3f", show[[v]]))
for (v in c("mean_S", "meanv2", "delta_t2", "delta_52", "true_sd", "est_sd")) show[[v]] <- sprintf("%.2f", show[[v]])
print(show, row.names = FALSE)
cat("\nmeanv2 = ||mean over reps of v||^2 (the real non-centrality, needs R large); delta_t2 / delta_52 = mean ||v(y*)||^2 with y* = pi(beta_tilde) / de-noised pi0-hat.",
    "\ntrue_sd = sd(eta_0); est_sd = R5's estimated signal sd.  wR* = sum w_ref / sum w(beta_hat).  Right reference: defl 1.00, size .05 everywhere.\n")
