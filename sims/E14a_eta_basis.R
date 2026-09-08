# E14a -- DOES THE BASIS BELONG ON THE PROBABILITY SCALE OR THE INDEX SCALE?
#
# E13 found the one failure region: at kappa = .25 the smooth (EDGE) statistic rejects .574 / .884
# of CORRECT models at prevalence .15 / .08, worse than the decile statistic, while at fixed p both
# are fine.  Recorded mechanism: the fitted probabilities crowd towards zero, so the EDGE polynomials
# -- built on the group-mean PROBABILITY pbar -- are compressed and mis-shaped.
#
# THE HYPOTHESIS TO TEST.  Shifting the intercept to lower the prevalence is, on the INDEX scale, a
# pure translation of eta.  poly(etabar + c, k) spans the same subspace as poly(etabar, k), so a basis
# built on the group-mean index is translation-EQUIVARIANT: the statistic should not see the shift at
# all.  On the probability scale the same shift is expit(), which is nonlinear, so poly(pbar, k)
# genuinely changes shape -- and the selection non-centrality, which E9/E9b found is almost entirely
# CUBIC in the group index, can rotate into the degree-2 basis that the kappa-regime rule uses.
# That would explain why the SMOOTH statistic (degree 2, g = .11 kappa) fails harder than the decile
# statistic (full space, g = .24 kappa): the aggregate inflation absorbs the term, the narrow one does not.
#
# WHAT IS MEASURED.  Same seed and same draw() as E13, so the data are IDENTICAL and the comparison is
# paired; the hl / edge columns must reproduce E13 exactly (a correctness check on this script).
#   (1) size of the decile statistic, of EDGE on pbar (the shipped statistic), and of EDGE on etabar;
#   (2) the decomposition that says WHY.  For each basis, mean(S) / mean(tr(P_Z Sigma_0)) is the ratio
#       the reference has to reproduce -- it is 1 when the reference is right and the inflation g is
#       exactly the missing part.  Splitting mean(S) into the selection mean m = E(v) projected into
#       the basis, ||P_Z m||^2, and the remainder separates a REFERENCE error (central part off 1)
#       from a SELECTION error (non-centrality share large).  The current inflation is g = .11 kappa
#       = .0275 here, so a non-centrality share far above that is the diagnosis.
# Usage: Rscript E14a_eta_basis.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500
W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260911                      # E13's seed -- do not change, the paired comparison depends on it
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))

.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-12, 12))$root
mk <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, G = 10, prev = 0.5) {
  Sig <- if (p == 5) diag(5) else ar^abs(outer(1:p, 1:p, "-"))
  b0 <- if (p == 5) c(.5, -.4, .3, -.3, .2) else { b <- rep_len(c(.35, -.3, .25, -.2, .15), p); b * norm / sqrt(drop(t(b) %*% Sig %*% b)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), G = G, prev = prev)
}
CELLS <- list(
  mk("k25_prev50", 400, 100, 416, prev = .50),   # the balanced control -- the eta basis must not break it
  mk("k25_prev30", 400, 100, 416, prev = .30),
  mk("k25_prev15", 400, 100, 416, prev = .15),   # E13: edge .574
  mk("k25_prev08", 400, 100, 416, prev = .08),   # E13: edge .884
  mk("k10_prev08", 400,  40, 416, prev = .08),   # is the failure kappa-graded?
  mk("A_prev50",   500,   5, 100, prev = .50),   # fixed-p controls, degree 3
  mk("A_prev08",   500,   5, 100, prev = .08))

draw <- function(d, r) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else .a_for(eta, d$prev)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
one <- function(r, d) {
  z <- draw(d, r); if (sum(z$y) < 2 * d$G || sum(1 - z$y) < 2 * d$G) return(NULL)
  pc <- tryCatch(cf_pieces(z$X, z$y, d$lambda, G = d$G), error = function(e) NULL); if (is.null(pc)) return(NULL)
  Sig0 <- tryCatch(cf_reference_R5(z$X, z$y, d$lambda, pc)$R5, error = function(e) NULL); if (is.null(Sig0)) return(NULL)
  kappa <- d$p / d$n
  ad <- cf_adaptive_degree(z$X, z$y, pc)
  kdeg <- if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L          # the shipped rho2k rule
  g <- kappa * c(.11, .11, 1.15)[kdeg]                            # basis-specific inflation, as shipped
  idx <- split(seq_len(d$n), pc$groups)
  pbar <- vapply(idx, function(I) mean(pc$pi[I]), 0.0)
  eta  <- drop(cbind(1, z$X) %*% pc$beta)
  ebar <- vapply(idx, function(I) mean(eta[I]), 0.0)
  Zp <- tryCatch(as.matrix(stats::poly(pbar, kdeg)), error = function(e) NULL)
  Ze <- tryCatch(as.matrix(stats::poly(ebar, kdeg)), error = function(e) NULL)
  if (is.null(Zp) || is.null(Ze)) return(NULL)
  tr_of <- function(Z) sum(diag(Z %*% solve(crossprod(Z), t(Z)) %*% Sig0))   # tr(P_Z Sigma_0)
  data.frame(rep = r, cell = d$name, n = d$n, p = d$p, prev = mean(z$y), deg = kdeg, rho_hat = ad$rho_hat,
             hl   = cf_pvalue(pc$S_dec, (1 + kappa * .24) * Sig0),
             edge = cf_pvalue(cf_qf(pc$v, Zp), (1 + g) * Sig0, Zp),      # shipped: basis on pbar
             eta  = cf_pvalue(cf_qf(pc$v, Ze), (1 + g) * Sig0, Ze),      # proposed: basis on etabar
             S_p = cf_qf(pc$v, Zp), S_e = cf_qf(pc$v, Ze), tr_p = tr_of(Zp), tr_e = tr_of(Ze),
             pbar_min = min(pbar), pbar_max = max(pbar), ebar_sd = sd(ebar),
             t(setNames(pc$v, paste0("v", seq_along(pc$v)))),
             t(setNames(pbar, paste0("pb", seq_along(pbar)))),
             t(setNames(ebar, paste0("eb", seq_along(ebar)))))
}
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c(ls(), ".a_for"))                       # ls() omits dot-prefixed names
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d), error = function(e) NULL)))
  cat(sprintf("done %-12s reps=%s (%.1f min)\n", d$name, if (is.null(o)) "0" else nrow(o),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl)
write.csv(res, file.path(OUT, "E14a_eta_basis_perrep.csv"), row.names = FALSE)

# ---- size, and the decomposition that explains it -----------------------
rej <- function(v) mean(v < .05, na.rm = TRUE)
sm <- do.call(rbind, lapply(CELLS, function(d) {
  A <- res[res$cell == d$name, , drop = FALSE]; if (!nrow(A)) return(NULL)
  V <- as.matrix(A[, paste0("v", seq_len(d$G))]); m <- colMeans(V)          # the selection mean E(v)
  nc <- function(pre) mean(vapply(seq_len(nrow(A)), function(i) {
    x <- as.numeric(A[i, paste0(pre, seq_len(d$G))]); Z <- as.matrix(stats::poly(x, A$deg[i]))
    drop(crossprod(crossprod(Z, m), solve(crossprod(Z), crossprod(Z, m)))) }, 0.0))
  data.frame(cell = d$name, kappa = d$p / d$n, prev = round(mean(A$prev), 3), deg = A$deg[1],
             rho_hat = round(mean(A$rho_hat), 2),
             size_hl = rej(A$hl), size_edge = rej(A$edge), size_eta = rej(A$eta),
             # ratio the reference must reproduce; 1 = correct.  split into central + selection parts
             ratio_p = mean(A$S_p) / mean(A$tr_p), ncsh_p = nc("pb") / mean(A$tr_p),
             ratio_e = mean(A$S_e) / mean(A$tr_e), ncsh_e = nc("eb") / mean(A$tr_e),
             g_used = round((d$p / d$n) * c(.11, .11, 1.15)[A$deg[1]], 3),
             pbar_rng = sprintf("%.3f-%.3f", mean(A$pbar_min), mean(A$pbar_max)), reps = nrow(A))
}))
write.csv(sm, file.path(OUT, "E14a_eta_basis_summary.csv"), row.names = FALSE)
options(width = 220)
cat(sprintf("\n=== E14a | probability-scale vs index-scale basis | R = %d | %.1f min ===\n",
            R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("size_edge = basis on group-mean PROBABILITY (shipped) | size_eta = basis on group-mean INDEX (proposed)\n")
cat("ratio = mean(S)/mean(tr(P_Z Sigma_0)); the inflation g_used is meant to be ratio - 1.  ncsh = the share of that ratio coming from the selection mean.\n\n")
print(sm[, c("cell", "kappa", "prev", "deg", "rho_hat", "size_hl", "size_edge", "size_eta", "ratio_p", "ncsh_p", "ratio_e", "ncsh_e", "g_used", "pbar_rng", "reps")],
      digits = 3, row.names = FALSE)
