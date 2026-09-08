# E14c -- IS THE BOUNDARY A CHECKABLE PRECONDITION RATHER THAN A PROHIBITION?
#
# E14a settled the mechanism.  The de-noised reference is right at every prevalence (the central part
# of mean(S)/mean(tr) sits at .88-1.00 from prevalence .50 down to .08), the basis scale is irrelevant,
# and the whole failure is the selection non-centrality: the O(kappa) term that the calibrated
# inflation absorbs when outcomes are balanced grows to nine times the reference trace at prevalence
# .08.  E14b confirms it by removing the selection entirely with an X-only grouping.
#
# That reframes the boundary.  The selection displacement enters the residual through the
# standardisation r_g = (sum y - sum pi) / sqrt(V_g), and V_g = sum pi(1-pi) over the group is what
# collapses when a group holds few events -- not the prevalence as such.  Prevalence, n, and G all act
# on the test only through that one quantity.  If the non-centrality share is a function of the
# WEAKEST group's V_g, then the user does not need a rule about prevalence at all: they need to check
# a number they can compute from their own fit before running the test, and can restore validity by
# choosing fewer groups.  Fewer groups is also what the clinical literature already does with few events.
#
# So: G in {5, 10, 20} crossed with prevalence in {.50, .30, .15, .08}, at kappa = .25 and at fixed p,
# recording min_g V_g and the expected events in the smallest group alongside the size and the
# non-centrality share.  Same seed and draw as E13/E14a/E14b.
# Usage: Rscript E14c_events_per_group.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500
W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260911
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))

.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-12, 12))$root
mk <- function(tag, n, p, lambda, G, prev) {
  Sig <- if (p == 5) diag(5) else 0.7^abs(outer(1:p, 1:p, "-"))
  b0 <- if (p == 5) c(.5, -.4, .3, -.3, .2) else { b <- rep_len(c(.35, -.3, .25, -.2, .15), p); b * 2.31 / sqrt(drop(t(b) %*% Sig %*% b)) }
  list(name = sprintf("%s_prev%02d_G%d", tag, round(100 * prev), G), tag = tag,
       n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), G = G, prev = prev)
}
CELLS <- c(
  unlist(lapply(c(.50, .30, .15, .08), function(pv) lapply(c(5, 10, 20), function(G) mk("k25", 400, 100, 416, G, pv))), recursive = FALSE),
  unlist(lapply(c(.50, .15, .08),      function(pv) lapply(c(5, 10, 20), function(G) mk("A",   500,   5, 100, G, pv))), recursive = FALSE))

draw <- function(d, r) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0)
  a <- if (abs(d$prev - .5) < 1e-8) 0 else .a_for(eta, d$prev)
  list(X = X, y = rbinom(d$n, 1, cf_expit(a + eta)))
}
GMAX <- 20                              # v / pbar are padded to GMAX so cells with different G can be rbound
pad <- function(x, pre) setNames(as.data.frame(t(c(x, rep(NA_real_, GMAX - length(x))))), paste0(pre, seq_len(GMAX)))
one <- function(r, d) {
  # guard at G, not 2G: at prevalence .08 with G = 20 the design is genuinely marginal, and dropping
  # those replications would hide the very failure mode the run is here to measure
  z <- draw(d, r); if (sum(z$y) < d$G || sum(1 - z$y) < d$G) return(NULL)
  pc <- tryCatch(cf_pieces(z$X, z$y, d$lambda, G = d$G), error = function(e) NULL); if (is.null(pc)) return(NULL)
  Sig0 <- tryCatch(cf_reference_R5(z$X, z$y, d$lambda, pc)$R5, error = function(e) NULL); if (is.null(Sig0)) return(NULL)
  kappa <- d$p / d$n
  ad <- cf_adaptive_degree(z$X, z$y, pc)
  kdeg <- if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L
  g <- kappa * c(.11, .11, 1.15)[kdeg]
  idx <- split(seq_len(d$n), pc$groups)
  w <- pmax(pc$pi * (1 - pc$pi), 1e-8)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)                     # the standardising weight per group
  ev <- vapply(idx, function(I) sum(pc$pi[I]), 0.0)                 # expected events per group
  pbar <- vapply(idx, function(I) mean(pc$pi[I]), 0.0)
  Zp <- tryCatch(as.matrix(stats::poly(pbar, kdeg)), error = function(e) NULL); if (is.null(Zp)) return(NULL)
  data.frame(rep = r, cell = d$name, tag = d$tag, G = d$G, prev = mean(z$y), deg = kdeg,
             hl   = cf_pvalue(pc$S_dec, (1 + kappa * .24) * Sig0),
             edge = cf_pvalue(cf_qf(pc$v, Zp), (1 + g) * Sig0, Zp),
             S_p = cf_qf(pc$v, Zp), tr_p = sum(diag(Zp %*% solve(crossprod(Zp), t(Zp)) %*% Sig0)),
             S_d = pc$S_dec, tr_d = sum(diag(Sig0)),
             Vmin = min(Vg), evmin = min(ev),
             pad(pc$v, "v"), pad(pbar, "pb"))
}
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c(ls(), ".a_for"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d), error = function(e) NULL)))
  cat(sprintf("done %-16s reps=%s (%.1f min)\n", d$name, if (is.null(o)) "0" else nrow(o),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
stopCluster(cl)
write.csv(res, file.path(OUT, "E14c_events_perrep.csv"), row.names = FALSE)

rej <- function(v) mean(v < .05, na.rm = TRUE)
sm <- do.call(rbind, lapply(CELLS, function(d) {
  A <- res[res$cell == d$name, , drop = FALSE]; if (!nrow(A)) return(NULL)
  m <- colMeans(as.matrix(A[, paste0("v", seq_len(d$G))]))
  ncp <- mean(vapply(seq_len(nrow(A)), function(i) {
    Z <- as.matrix(stats::poly(as.numeric(A[i, paste0("pb", seq_len(d$G))]), A$deg[i]))
    drop(crossprod(crossprod(Z, m), solve(crossprod(Z), crossprod(Z, m)))) }, 0.0))
  data.frame(cell = d$name, tag = d$tag, G = d$G, prev = round(mean(A$prev), 3),
             evmin = round(mean(A$evmin), 1), Vmin = round(mean(A$Vmin), 2),
             size_hl = rej(A$hl), size_edge = rej(A$edge),
             ratio_d = mean(A$S_d) / mean(A$tr_d), ratio_e = mean(A$S_p) / mean(A$tr_p),
             ncsh_e = ncp / mean(A$tr_p), reps = nrow(A))
}))
sm <- sm[order(sm$tag, sm$G, -sm$prev), ]
write.csv(sm, file.path(OUT, "E14c_events_summary.csv"), row.names = FALSE)
options(width = 220)
cat(sprintf("\n=== E14c | is the boundary a precondition on the weakest group? | R = %d | %.1f min ===\n",
            R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("evmin = expected events in the smallest group; Vmin = its standardising weight sum.\n")
cat("If the size and the non-centrality share are governed by evmin/Vmin rather than by prevalence, the boundary is checkable and G is the dial.\n\n")
print(sm, digits = 3, row.names = FALSE)
cat("\n--- ordered by evmin, both designs pooled: does the failure track the weakest group? ---\n")
print(sm[order(sm$evmin), c("cell", "G", "prev", "evmin", "Vmin", "size_hl", "size_edge", "ncsh_e")], digits = 3, row.names = FALSE)
