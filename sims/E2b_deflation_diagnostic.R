# =====================================================================
# E2b -- WHERE does the 0.80 deflation come from, and WHAT observable tracks it?
#
# Under H0^str, per replication, compute the corrected statistics with three
# groupings held on the SAME fit:
#   (i)   random : deciles of pi_hat            (the test as used)
#   (ii)  oracle : deciles of the TRUE index eta_0 = X beta_0   (fixed, not data-driven)
#   (iii) debias : deciles of the debiased index X beta_tilde
# and record E[S_c] / tr(Omega_MLE) for each -- the deflation ratio.  If the
# ratio is ~1 under (ii), the deflation is a grouping effect; if it stays ~0.8,
# it lives in mu_hat / the residual geometry and the grouping is innocent.
#
# Alongside, the Bellec (2022) observables from the engine, so that candidate
# closed-form deflation factors can be compared with the measured ratio.
#
# Cells: the six registered ones (A: lambda 100/137/200; B: lambda 50/416/1000)
# plus a kappa sweep at lambda = 416 with the design-B recipe: p in {40, 100, 160}
# (kappa = 0.10, 0.25, 0.40), beta cycling the same five values scaled to ||beta|| = 2.31.
#
# Usage: Rscript E2b_deflation_diagnostic.R [R] [workers]
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
  else {
    beta0 <- rep_len(c(0.35, -0.30, 0.25, -0.20, 0.15), p); beta0 <- beta0 * beta_norm / sqrt(sum(beta0^2))
    Sc <- chol(rho^abs(outer(1:p, 1:p, "-")))
  }
  list(tag = tag, n = n, p = p, lambda = lambda, beta0 = beta0, Sc = Sc)
}
CELLS <- list(
  design("A", 500,   5, 100,  0,   NA), design("A", 500,   5, 137,  0,   NA), design("A", 500,   5, 200,  0,   NA),
  design("B", 400, 100,  50,  0.7, 2.31), design("B", 400, 100, 416, 0.7, 2.31), design("B", 400, 100, 1000, 0.7, 2.31),
  design("B", 400,  40, 416,  0.7, 2.31), design("B", 400, 160, 416, 0.7, 2.31))
names(CELLS) <- vapply(CELLS, function(d) sprintf("%s_n%d_p%d_lam%g", d$tag, d$n, d$p, d$lambda), "")

grp <- function(v, G = 10) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)
edge_df <- function(Om, Z) { e <- eigen(Om, symmetric = TRUE); Sh <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors)); P <- Z %*% solve(crossprod(Z), t(Z)); sum(diag(Sh %*% P %*% Sh)) }

one_rep <- function(r, d) {
  set.seed(SEED + r)
  X <- matrix(rnorm(d$n * d$p), d$n, d$p) %*% d$Sc
  eta0 <- drop(X %*% d$beta0); y <- rbinom(d$n, 1, cf_expit(eta0))
  X1 <- cbind(1, X)
  beta <- cf_ridge(X1, y, d$lambda, diag(c(0, rep(1, d$p))))
  pc_r <- cf_pieces(X, y, d$lambda, beta = beta)                              # (i) random
  pc_o <- cf_pieces(X, y, d$lambda, beta = beta, groups = grp(eta0))          # (ii) oracle
  pc_d <- cf_pieces(X, y, d$lambda, beta = beta, groups = grp(drop(X1 %*% pc_r$beta_tilde)))  # (iii) debiased
  row <- function(pc, lab) data.frame(rep = r, cell = names(CELLS)[sapply(CELLS, identical, d)], grouping = lab,
    S_dec = pc$S_dec, S_edge = pc$S_edge, tr_OmMLE = pc$obs["tr_OmMLE"], df_edge = if (is.null(pc$Z)) NA else edge_df(pc$Om_MLE, pc$Z),
    S_dec_unc = pc$S_dec_unc, tr_OmK = pc$obs["tr_OmK"],
    p_dec = cf_pvalue(pc$S_dec, pc$Om_MLE), p_edge = cf_pvalue(pc$S_edge, pc$Om_MLE, pc$Z),
    t(pc$obs[c("df_M", "trV_M", "v_hat", "r2_hat", "gamma_hat", "df_F", "trV_F", "trW", "norm_bhat", "norm_btilde")]))
  rbind(row(pc_r, "random"), row(pc_o, "oracle"), row(pc_d, "debiased"))
}

t0 <- Sys.time()
cl <- makeCluster(WORKERS); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one_rep, d = d))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E2b_deflation_perrep.csv"), row.names = FALSE)

# ---- summary --------------------------------------------------------------
agg <- aggregate(cbind(S_dec, S_edge, tr_OmMLE, df_edge, tr_OmK, S_dec_unc, df_M, trV_M, v_hat, r2_hat, gamma_hat, df_F, trV_F, trW, norm_bhat, norm_btilde,
                       rej_dec = p_dec < 0.05, rej_edge = p_edge < 0.05) ~ cell + grouping, data = res, FUN = mean)
agg$defl_dec  <- agg$S_dec  / agg$tr_OmMLE
agg$defl_edge <- agg$S_edge / agg$df_edge
# candidate observable ratios (computed on the random-grouping rows; geometry-only, same across groupings)
agg$n <- as.integer(sub(".*_n([0-9]+)_.*", "\1", agg$cell)); agg$p <- as.integer(sub(".*_p([0-9]+)_.*", "\1", agg$cell))
# candidate OBSERVABLE deflation factors (< 1 when the ridge fit spends fewer df than the MLE)
agg$cand_trV <- agg$trV_F / agg$trV_M                       # tr V_F / tr V_M
agg$cand_df  <- (agg$n - agg$df_F) / (agg$n - agg$df_M)     # residual df of the MLE fit / of the ridge fit
agg$cand_r2  <- agg$r2_hat / (agg$trW / agg$n)              # observed residual second moment / naive noise level
agg$R <- R
agg <- agg[order(agg$cell, agg$grouping), ]
write.csv(agg, file.path(OUT, "E2b_deflation_summary.csv"), row.names = FALSE)

options(width = 220)
cat(sprintf("\n=== E2b | R = %d per cell | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
show <- agg[, c("cell", "grouping", "defl_dec", "defl_edge", "rej_dec", "rej_edge", "tr_OmMLE", "df_M", "df_F", "cand_trV", "cand_df", "cand_r2", "norm_btilde")]
for (v in c("defl_dec", "defl_edge", "cand_trV", "cand_df", "cand_r2")) show[[v]] <- sprintf("%.3f", show[[v]])
for (v in c("rej_dec", "rej_edge")) show[[v]] <- sprintf("%.3f", show[[v]])
for (v in c("tr_OmMLE", "df_M", "norm_btilde")) show[[v]] <- sprintf("%.2f", show[[v]])
print(show, row.names = FALSE)
cat("\ndefl_* = mean S_c / mean tr(reference covariance) -- 1.00 means first-order theory is exact; rej_* = analytic (deflate = 1) size at 5%.",
    "\ncand_trV = trV_F/trV_M, cand_df = (n-df_F)/(n-df_M), cand_r2 = r2_hat/(trW/n): candidate OBSERVABLE deflation factors, to compare with defl_*.\n")
