# E14d -- WHAT DOES THE y-FREE GROUPING COST?
#
# E14b showed that grouping on a function of the predictors alone removes the low-prevalence failure
# outright: at kappa = .25, prevalence .08, the smooth statistic goes from rejecting .884 of correct
# models to .066, and the ratio the reference must reproduce returns to 1.00 / 1.18.  That proves the
# de-noised reference was never the problem.  It does NOT yet make X-only grouping a proposal, because
# the whole point of grouping on the fitted index is that the groups are ordered by the quantity whose
# calibration is being tested; grouping on the first principal component of X orders them by something
# else, and a departure that bends the calibration curve along the fitted index may average out.
#
# This run decides whether that worry is real, by putting the two groupings on the SAME data under the
# SAME departures used in Table 5 (E4d): a quadratic bend in the index, and an off-index interaction.
# Prevalence is held at the target AFTER the departure is added, so the departure and the prevalence
# do not move together.  At prevalence .08 the fitted grouping is invalid (size .88), so its power
# column there is not a power at all and is reported only to make that visible.
# Usage: Rscript E14d_power_under_yfree_grouping.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 400
W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260912
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))

.a_for <- function(eta, target) stats::uniroot(function(a) mean(cf_expit(a + eta)) - target, c(-15, 15))$root
mk <- function(name, n, p, lambda, G = 10, prev = 0.5) {
  Sig <- if (p == 5) diag(5) else 0.7^abs(outer(1:p, 1:p, "-"))
  b0 <- if (p == 5) c(.5, -.4, .3, -.3, .2) else { b <- rep_len(c(.35, -.3, .25, -.2, .15), p); b * 2.31 / sqrt(drop(t(b) %*% Sig %*% b)) }
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = chol(Sig), G = G, prev = prev,
       sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (p == 5) 1 else sqrt(1 + .49))
}
CELLS <- list(mk("k25_prev50", 400, 100, 416, prev = .50),    # both groupings valid here -- the fair comparison
              mk("k25_prev08", 400, 100, 416, prev = .08))    # fitted is invalid here; shown to make that visible
DEPS <- rbind(data.frame(dep = "null", rel = 0),
              data.frame(dep = "quadratic", rel = c(.5, 1, 1.5)),
              data.frame(dep = "offindex", rel = c(1, 1.5)))
grp <- function(v, G) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)

one <- function(r, d, dep, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc
  e0 <- drop(X %*% d$b0); Zs <- e0 / d$sd_eta
  e0 <- e0 + switch(dep, null = 0,
                    quadratic = d$sd_eta * rel * (Zs^2 - 1) / sqrt(2),
                    offindex  = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  a <- .a_for(e0, d$prev)                        # prevalence held AFTER the departure
  y <- rbinom(d$n, 1, cf_expit(a + e0))
  if (sum(y) < d$G || sum(1 - y) < d$G) return(NULL)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, d$p))); kappa <- d$p / d$n
  gset <- tryCatch(list(fitted = grp(drop(X1 %*% cf_ridge(X1, y, d$lambda, Dpen)), d$G),
                        pc1    = grp(drop(X %*% eigen(crossprod(scale(X, scale = FALSE)) / d$n, symmetric = TRUE)$vectors[, 1]), d$G)),
                   error = function(e) NULL)
  if (is.null(gset)) return(NULL)
  do.call(rbind, lapply(names(gset), function(gn) {
    pc <- tryCatch(cf_pieces(X, y, d$lambda, G = d$G, groups = gset[[gn]]), error = function(e) NULL); if (is.null(pc)) return(NULL)
    Sig0 <- tryCatch(cf_reference_R5(X, y, d$lambda, pc)$R5, error = function(e) NULL); if (is.null(Sig0)) return(NULL)
    ad <- cf_adaptive_degree(X, y, pc); kdeg <- if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L
    g <- kappa * c(.11, .11, 1.15)[kdeg]
    idx <- split(seq_len(d$n), pc$groups)
    pbar <- vapply(idx, function(I) mean(pc$pi[I]), 0.0)
    Zp <- tryCatch(as.matrix(stats::poly(pbar, kdeg)), error = function(e) NULL); if (is.null(Zp)) return(NULL)
    data.frame(rep = r, cell = d$name, dep = dep, rel = rel, grouping = gn, prev = mean(y),
               hl = cf_pvalue(pc$S_dec, (1 + kappa * .24) * Sig0),
               edge = cf_pvalue(cf_qf(pc$v, Zp), (1 + g) * Sig0, Zp))
  }))
}
t0 <- Sys.time(); cl <- makeCluster(W)
clusterExport(cl, c(ls(), ".a_for"))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(j) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d, DEPS$dep[j], DEPS$rel[j]), error = function(e) NULL)))
  cat(sprintf("done %-12s %-10s rel=%.1f rows=%s (%.1f min)\n", d$name, DEPS$dep[j], DEPS$rel[j],
              if (is.null(o)) "0" else nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))))
stopCluster(cl)
write.csv(res, file.path(OUT, "E14d_power_yfree_perrep.csv"), row.names = FALSE)

rej <- function(v) mean(v < .05, na.rm = TRUE)
sm <- aggregate(cbind(hl, edge) ~ cell + dep + rel + grouping, res, rej)
sm <- sm[order(sm$cell, sm$dep, sm$rel, sm$grouping), ]
write.csv(sm, file.path(OUT, "E14d_power_yfree_summary.csv"), row.names = FALSE)
options(width = 200)
cat(sprintf("\n=== E14d | what the X-only grouping costs in power | R = %d | %.1f min ===\n",
            R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("At prevalence .50 both groupings hold their level, so the comparison is a fair one.\n")
cat("At prevalence .08 the fitted grouping rejects ~.88 of CORRECT models, so its numbers there are not powers.\n\n")
for (cc in unique(sm$cell)) {
  cat("---", cc, "---\n")
  A <- sm[sm$cell == cc, ]
  w <- reshape(A[, c("dep", "rel", "grouping", "hl", "edge")], idvar = c("dep", "rel"),
               timevar = "grouping", direction = "wide")
  print(w, digits = 3, row.names = FALSE); cat("\n")
}
