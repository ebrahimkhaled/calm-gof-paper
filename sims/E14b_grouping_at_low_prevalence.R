# E14b -- IF THE SELECTION MEAN IS THE WHOLE STORY, A y-FREE GROUPING SHOULD FIX IT.
#
# E14a killed the basis hypothesis and replaced it with a sharper one.  Decomposing the ratio that
# the reference has to reproduce, mean(S) / mean(tr(P_Z Sigma_0)), into the part carried by the
# selection mean m = E(v) and the rest, the central part sits near 1 at every prevalence -- the
# de-noised reference is RIGHT -- while the non-centrality share runs away as the outcome unbalances.
# The basis scale is irrelevant to that: putting the polynomials on the group-mean index instead of
# the group-mean probability moves the ratio the wrong way.
#
# So the failure is not the basis and not the variance.  It is that the groups are ranked on an index
# that depends on y, and at low prevalence that dependence is violently amplified: a unit in a bottom
# group with y = 1 raises its own pi_hat and leaves the group, so the bottom groups run systematically
# short of events, and the standardisation by sqrt(V_g) -- with V_g = sum pi(1-pi) tiny when pi ~ .02 --
# divides that displacement by a small number.
#
# E5 already built the groupings that test this, and the engine already exposes two of them:
#   fitted   deciles of eta_hat(lambda)          -- the shipped default, self-influence O(kappa)
#   heavyx5  deciles of eta_hat(5 lambda)        -- grouping = "heavy" in calm.gof(), E5's fix
#   heavyx20 deciles of eta_hat(20 lambda)       -- more of the same
#   pc1      deciles of the first PC of X        -- a function of X ALONE, so zero selection by
#            construction.  Not a proposal (it groups on the wrong thing and will cost power); it is
#            the diagnostic that proves whether selection is the entire mechanism.
# Same seed and draw as E13/E14a, so every comparison is paired on identical data.
# Usage: Rscript E14b_grouping_at_low_prevalence.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500
W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260911
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
  mk("k25_prev50", 400, 100, 416, prev = .50),
  mk("k25_prev15", 400, 100, 416, prev = .15),
  mk("k25_prev08", 400, 100, 416, prev = .08),
  mk("k10_prev08", 400,  40, 416, prev = .08),
  mk("A_prev08",   500,   5, 100, prev = .08))
GRPS <- c("fitted", "heavyx5", "heavyx20", "pc1")

draw <- function(d, r) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else .a_for(eta, d$prev)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
grp <- function(v, G) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)
one <- function(r, d) {
  z <- draw(d, r); if (sum(z$y) < 2 * d$G || sum(1 - z$y) < 2 * d$G) return(NULL)
  X1 <- cbind(1, z$X); Dpen <- diag(c(0, rep(1, d$p))); kappa <- d$p / d$n
  gset <- tryCatch(list(
    fitted   = grp(drop(X1 %*% cf_ridge(X1, z$y,      d$lambda, Dpen)), d$G),
    heavyx5  = grp(drop(X1 %*% cf_ridge(X1, z$y,  5 * d$lambda, Dpen)), d$G),
    heavyx20 = grp(drop(X1 %*% cf_ridge(X1, z$y, 20 * d$lambda, Dpen)), d$G),
    pc1      = grp(drop(z$X %*% eigen(crossprod(scale(z$X, scale = FALSE)) / d$n, symmetric = TRUE)$vectors[, 1]), d$G)),
    error = function(e) NULL)
  if (is.null(gset)) return(NULL)
  do.call(rbind, lapply(names(gset), function(gn) {
    pc <- tryCatch(cf_pieces(z$X, z$y, d$lambda, G = d$G, groups = gset[[gn]]), error = function(e) NULL); if (is.null(pc)) return(NULL)
    Sig0 <- tryCatch(cf_reference_R5(z$X, z$y, d$lambda, pc)$R5, error = function(e) NULL); if (is.null(Sig0)) return(NULL)
    ad <- cf_adaptive_degree(z$X, z$y, pc)
    kdeg <- if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L
    g <- kappa * c(.11, .11, 1.15)[kdeg]
    idx <- split(seq_len(d$n), pc$groups)
    pbar <- vapply(idx, function(I) mean(pc$pi[I]), 0.0)
    Zp <- tryCatch(as.matrix(stats::poly(pbar, kdeg)), error = function(e) NULL); if (is.null(Zp)) return(NULL)
    data.frame(rep = r, cell = d$name, grouping = gn, prev = mean(z$y), deg = kdeg,
               hl   = cf_pvalue(pc$S_dec, (1 + kappa * .24) * Sig0),
               edge = cf_pvalue(cf_qf(pc$v, Zp), (1 + g) * Sig0, Zp),
               S_p = cf_qf(pc$v, Zp), tr_p = sum(diag(Zp %*% solve(crossprod(Zp), t(Zp)) %*% Sig0)),
               S_d = pc$S_dec, tr_d = sum(diag(Sig0)),
               t(setNames(pc$v, paste0("v", seq_along(pc$v)))),
               t(setNames(pbar, paste0("pb", seq_along(pbar)))))
  }))
}
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c(ls(), ".a_for"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d), error = function(e) NULL)))
  cat(sprintf("done %-12s rows=%s (%.1f min)\n", d$name, if (is.null(o)) "0" else nrow(o),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl)
write.csv(res, file.path(OUT, "E14b_grouping_perrep.csv"), row.names = FALSE)

rej <- function(v) mean(v < .05, na.rm = TRUE)
sm <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(GRPS, function(gn) {
  A <- res[res$cell == d$name & res$grouping == gn, , drop = FALSE]; if (!nrow(A)) return(NULL)
  m <- colMeans(as.matrix(A[, paste0("v", seq_len(d$G))]))            # the selection mean E(v)
  ncp <- mean(vapply(seq_len(nrow(A)), function(i) {
    Z <- as.matrix(stats::poly(as.numeric(A[i, paste0("pb", seq_len(d$G))]), A$deg[i]))
    drop(crossprod(crossprod(Z, m), solve(crossprod(Z), crossprod(Z, m)))) }, 0.0))
  data.frame(cell = d$name, kappa = d$p / d$n, prev = round(mean(A$prev), 3), grouping = gn,
             size_hl = rej(A$hl), size_edge = rej(A$edge),
             ratio_d = mean(A$S_d) / mean(A$tr_d), ratio_e = mean(A$S_p) / mean(A$tr_p),
             ncsh_e = ncp / mean(A$tr_p), norm_m = sqrt(sum(m^2)), reps = nrow(A))
}))))
write.csv(sm, file.path(OUT, "E14b_grouping_summary.csv"), row.names = FALSE)
options(width = 220)
cat(sprintf("\n=== E14b | does a y-free grouping remove the low-prevalence failure? | R = %d | %.1f min ===\n",
            R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("ratio = mean(S)/mean(tr) under the UNinflated reference (1 = right).  ncsh_e = share of it carried by the selection mean.\n")
cat("norm_m = ||E(v)||, the size of the selection displacement itself.  pc1 has zero selection BY CONSTRUCTION -- it is the control, not a proposal.\n\n")
print(sm, digits = 3, row.names = FALSE)
