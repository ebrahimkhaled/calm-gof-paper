# E11 -- can the selection non-centrality be estimated PER DATASET instead of by the fitted constants 0.24*kappa / 1.15*kappa?
# Motivation (E10/E10b): ||E v||^2 is NOT a function of kappa. At kappa = .25 it is ~7x larger under AR 0.3 than AR 0.7 and ~2x larger at
# n = 800 than n = 400. The mechanism (Section 4.3) says the term is created by observation i's own y_i pushing eta_hat_i across a group
# boundary. That displacement is observable: delta_i = h_ii (y_i - pi_hat_i) / (w_i (1 - h_ii)) is the one-step leave-one-out correction
# (Pregibon 1981). So compute the group sums with and without it and take the difference:
#     D_g = V_gg^{-1/2} [ sum_{i: group(eta_hat)=g} (y_i - pi_hat_i)  -  sum_{i: group(eta_hat - delta)=g} (y_i - pi_hat_i) ]
# D is an observable, per-dataset estimate of E[v] (noisy). Two candidate references:
#   (a) inflation  g_obs = max(0, ||D||^2 - trace-based noise floor) / tr(Sigma5)     [same form as CALM, no fitted constant]
#   (b) non-central: refer S to the law of ||N(D, Sigma5)||^2 (Davies with non-centrality)  -- sharper if D is accurate
# Also record the two candidate SCALING LAWS for the fitted constant: kappa (current) and n*hbar^2/sd(eta_hat)^2 (the boundary-density
# argument: E[v_g] ~ sqrt(nG/wbar) * hbar * phi_g / sd(eta_hat)).
# Cells: the 4 calibration cells + the 4 held-out failures. Null only (this is about size). Usage: Rscript E11_observable_selection.R [R] [W]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args) >= 1) as.integer(args[1]) else 400; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260910
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(name, n, p, lambda, ar = 0.7, norm = 2.31) { Sig <- ar^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0)); list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(mk("A_k01", 500, 5, 100, ar = 0, norm = sqrt(sum(c(.5,-.4,.3,-.3,.2)^2))), mk("B_k10", 400, 40, 416), mk("B_k25", 400, 100, 416), mk("B_k40", 400, 160, 416),
              mk("H5_k25_ar03", 400, 100, 416, ar = 0.3), mk("H8_n800_k25", 800, 200, 832), mk("H3_k40_lam50", 400, 160, 50), mk("H1_k15", 400, 60, 416))
CELLS[[1]]$Sc <- diag(5); CELLS[[1]]$b0 <- c(.5, -.4, .3, -.3, .2)
G <- 10
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, cf_expit(drop(X %*% d$b0)))
  n <- d$n; X1 <- cbind(1, X); D <- diag(c(0, rep(1, d$p)))
  pc <- cf_pieces(X, y, d$lambda, G = G); R5 <- cf_reference_R5(X, y, d$lambda, pc)$R5
  bh <- pc$beta; eh <- drop(X1 %*% bh); ph <- cf_expit(eh); wh <- pmax(ph * (1 - ph), 1e-8); psi <- y - ph
  M <- crossprod(X1, wh * X1) + d$lambda * D
  h <- wh * rowSums((X1 %*% solve(M)) * X1)                       # self-influence h_ii
  delta <- h * psi / (wh * pmax(1 - h, 1e-6))                     # one-step LOO displacement of the index
  grp <- function(e) pmin(ceiling(rank(e, ties.method = "first") / (n / G)), G)
  g1 <- grp(eh); g0 <- grp(eh - delta)
  Vg <- as.numeric(tapply(wh, factor(g1, levels = 1:G), sum))
  s1 <- as.numeric(tapply(psi, factor(g1, levels = 1:G), sum)); s0 <- sapply(1:G, function(k) sum(psi[g0 == k]))
  Dvec <- (s1 - s0) / sqrt(Vg)                                     # observable estimate of E[v]
  kappa <- d$p / n; hbar <- mean(h); sde <- sd(eh)
  law_obs <- n * hbar^2 / sde^2                                    # boundary-density scaling law
  tr5 <- sum(diag(R5))
  g_obs <- sum(Dvec^2) / tr5
  # noise floor for ||D||^2: D is a difference of group sums over the few observations that changed group
  nch <- sum(g1 != g0); floor_est <- nch * mean(wh) / mean(Vg)
  data.frame(rep = r, cell = d$name, kappa = kappa, hbar = hbar, sd_eh = sde, law_obs = law_obs, tr5 = tr5,
             normD2 = sum(Dvec^2), g_obs = g_obs, nchanged = nch, floor_est = floor_est,
             S_dec = pc$S_dec, v2 = sum(pc$v^2),
             p_none = cf_pvalue(pc$S_dec, R5),
             p_kappa = cf_pvalue(pc$S_dec, (1 + 0.24 * kappa) * R5),
             p_obs = cf_pvalue(pc$S_dec, (1 + max(0, g_obs)) * R5),
             p_obs_sh = cf_pvalue(pc$S_dec, (1 + max(0, (sum(Dvec^2) - floor_est) / tr5)) * R5),
             p_law = cf_pvalue(pc$S_dec, (1 + 0.033 * law_obs) * R5),   # 0.033 chosen so that at B_k25 it reproduces 0.24*0.25 = .06
             vg = paste(round(pc$v, 3), collapse = "|"), Dg = paste(round(Dvec, 3), collapse = "|"))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d), error = function(e) NULL)))
  cat(sprintf("done %s n=%d (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl); write.csv(res[, setdiff(names(res), c("vg", "Dg"))], file.path(OUT, "E11_observable_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(kappa, hbar, sd_eh, law_obs, tr5, normD2, g_obs, nchanged, floor_est) ~ cell, res, mean)
sz <- aggregate(cbind(p_none, p_kappa, p_obs, p_obs_sh, p_law) ~ cell, res, rej)
# true ||E v||^2 (mean of v across reps) for comparison
tv <- do.call(rbind, lapply(split(res, res$cell), function(d) { V <- do.call(rbind, lapply(strsplit(d$vg, "\\|"), as.numeric))
  Dm <- do.call(rbind, lapply(strsplit(d$Dg, "\\|"), as.numeric)); data.frame(cell = d$cell[1], true_Ev2 = sum(colMeans(V)^2), meanD2 = sum(colMeans(Dm)^2), cor_shape = cor(colMeans(V), colMeans(Dm))) }))
agg <- merge(merge(agg, sz, by = "cell"), tv, by = "cell"); agg$g_true <- agg$true_Ev2 / agg$tr5; agg$R <- R
write.csv(agg, file.path(OUT, "E11_observable_summary.csv"), row.names = FALSE)
options(width = 250); cat(sprintf("\n=== E11 | observable selection non-centrality | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(agg[, c("cell", "kappa", "law_obs", "true_Ev2", "meanD2", "cor_shape", "g_true", "g_obs", "nchanged", "p_none", "p_kappa", "p_obs", "p_obs_sh", "p_law")], digits = 3, row.names = FALSE)
cat("\ng_true = ||E v||^2/tr(Sigma5) (the target); g_obs = mean per-dataset ||D||^2/tr; p_* = size at nominal .05 under each inflation.\n")
