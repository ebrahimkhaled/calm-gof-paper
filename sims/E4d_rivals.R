# E4d -- RIVAL tests on the SAME data as E4b/E4c/E4c2 (SEED + r regenerates X, eta, y exactly; merge by rep/cell/dep/rel).
# The author: "not me compare to me". Rivals a practitioner or a Series B referee would name:
#   grp        Janková, Shah, Bühlmann & Samworth (2020) GRP test (GRPtests::GRPtest, binomial, default 5 splits) -- the Series B rival;
#              tests a different null (model correct as a function of X) with its own lasso fit.
#   mle_hl     classical Hosmer-Lemeshow on the UNPENALISED maximum-likelihood fit, exact Moore-Spruill law (Omega_MLE), G = 10;
#   mle_edge   EDGE-3 on the MLE fit.  Both NA when the MLE does not exist (separation: |eta| > 30 or the IRLS diverges).
#   unc_hl     uncorrected HL on the RIDGE fit (the textbook test everyone runs) -- INVALID; reported as size + size-adjusted power.
#   cvhl       HL on 10-fold cross-validated ridge predictions (fixed lambda), chi-square(G-2) reference -- what careful practitioners do.
#   calslope   calibration intercept/slope LR test on the ridge fit (Cox 1958 recalibration; Van Calster 2016 "weak calibration"):
#              y ~ a + b*eta_hat, LR test of (a,b) = (0,1), chi-square(2).
#   spieg      Spiegelhalter (1986) z-test of calibration on the ridge fit.
#   uss        Copas (1989) unweighted sum of squares on the ridge fit with the Hosmer et al. (1997) moments (ML hat matrix).
# Cells: A, kappa .10/.25/.40; deps: null + quadratic .5/1/1.5 + offindex 1/1.5.   Usage: Rscript E4d_rivals.R [R] [workers]
suppressPackageStartupMessages({library(parallel); library(GRPtests)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 300; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260908
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416), mk("B", 400, 160, 416))
DEPS <- rbind(data.frame(dep = "null", rel = 0), data.frame(dep = "quadratic", rel = c(.5, 1, 1.5)), data.frame(dep = "offindex", rel = c(1, 1.5)))
G <- 10
hl_stat <- function(y, pi, G = 10) {            # textbook HL on given probabilities: sum (O-E)^2 / (E(1-pi_bar)), G equal-frequency groups
  g <- pmin(ceiling(rank(pi, ties.method = "first") / (length(y) / G)), G)
  O <- tapply(y, g, sum); E <- tapply(pi, g, sum); nb <- tapply(pi, g, length); pb <- E / nb
  sum((O - E)^2 / (nb * pb * (1 - pb)))
}
one <- function(r, d, dep, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, null = 0, quadratic = d$sd_eta * rel * (Zs^2 - 1) / sqrt(2), offindex = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta)); n <- d$n; X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, d$p)))
  # --- ridge fit (the deployed model) ---
  bh <- cf_ridge(X1, y, d$lambda, Dpen); eh <- drop(X1 %*% bh); ph <- cf_expit(eh); wh <- ph * (1 - ph)
  pc <- cf_pieces(X, y, d$lambda, G = G); unc_hl <- cf_pvalue(pc$S_dec_unc, pc$Om_MLE)
  # calibration intercept/slope LR test (chi2_2): free (a,b) vs offset(eta_hat)
  m1 <- suppressWarnings(glm(y ~ eh, family = binomial)); dev0 <- -2 * sum(y * eh - log1p(exp(eh)))
  calslope <- pchisq(dev0 - m1$deviance, 2, lower.tail = FALSE)
  # Spiegelhalter z
  dd <- 1 - 2 * ph; z <- sum((y - ph) * dd) / sqrt(sum(dd^2 * wh)); spieg <- 2 * pnorm(-abs(z))
  # Copas USS with Hosmer et al. (1997) moments, ML-form hat matrix at the ridge fit
  Fm <- crossprod(X1, wh * X1); dW <- dd * wh
  vU <- sum(dd^2 * wh) - drop(t(dW) %*% X1 %*% solve(Fm, crossprod(X1, dW)))
  uss <- 2 * pnorm(-abs((sum((y - ph)^2) - sum(wh)) / sqrt(max(vU, 1e-12))))
  # 10-fold CV ridge predictions, textbook HL chi2_{G-2}
  fold <- sample(rep_len(1:10, n)); pcv <- numeric(n)
  for (k in 1:10) { tr <- fold != k; bk <- cf_ridge(X1[tr, ], y[tr], d$lambda * mean(tr), Dpen); pcv[!tr] <- cf_expit(drop(X1[!tr, ] %*% bk)) }
  cvhl <- pchisq(hl_stat(y, pcv, G), G - 2, lower.tail = FALSE)
  # classical tests on the MLE (lambda = 0), NA when the MLE does not exist
  mle_hl <- mle_edge <- NA_real_; mle_ok <- FALSE
  # existence check via glm: converged, no fitted probability at the boundary, coefficients finite and moderate (Candès & Sur 2020)
  try({ gm <- suppressWarnings(glm.fit(X1, y, family = binomial())); pg <- gm$fitted.values
        if (gm$converged && all(pg > 1e-6 & pg < 1 - 1e-6) && max(abs(gm$coefficients)) < 20) {
          pm <- cf_pieces(X, y, 1e-8, G = G); mle_ok <- TRUE; mle_hl <- cf_pvalue(pm$S_dec_unc, pm$Om_MLE); mle_edge <- cf_pvalue(pm$S_edge_unc, pm$Om_MLE, pm$Z) } }, silent = TRUE)
  # GRP test (own lasso fit; nsplits default)
  grp <- tryCatch(as.numeric(GRPtest(X, y, fam = "binomial")), error = function(e) NA_real_)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, unc_hl = unc_hl, calslope = calslope, spieg = spieg, uss = uss, cvhl = cvhl,
             mle_ok = mle_ok, mle_hl = mle_hl, mle_edge = mle_edge, grp = grp)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(GRPtests)})))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) {
  out <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d, DEPS$dep[k], DEPS$rel[k]), error = function(e) NULL)))
  write.csv(out, file.path(OUT, sprintf("E4d_perrep_%s_%s_%g.csv", d$name, DEPS$dep[k], DEPS$rel[k])), row.names = FALSE)
  cat(sprintf("done %s %s rel %.1f  n=%d  (%.0f min)\n", d$name, DEPS$dep[k], DEPS$rel[k], nrow(out), as.numeric(difftime(Sys.time(), t0, units = "mins")))); out }))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E4d_rivals_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
vars <- c("unc_hl", "calslope", "spieg", "uss", "cvhl", "mle_hl", "mle_edge", "grp")
agg <- aggregate(res[, vars], by = res[, c("cell", "dep", "rel")], FUN = rej)
agg$mle_exists <- aggregate(mle_ok ~ cell + dep + rel, res, mean)$mle_ok
# size-adjusted power: reject at the empirical 5% quantile of the SAME test's null p-values in the SAME cell (for tests that do not hold the level)
adj <- do.call(rbind, lapply(split(res, res$cell), function(dc) { nul <- dc[dc$dep == "null", ]; do.call(rbind, lapply(split(dc, list(dc$dep, dc$rel), drop = TRUE), function(dk) {
  o <- data.frame(cell = dk$cell[1], dep = dk$dep[1], rel = dk$rel[1]); for (v in vars) { cv <- quantile(nul[[v]], .05, na.rm = TRUE); o[[paste0(v, "_adj")]] <- mean(dk[[v]] <= cv, na.rm = TRUE) }; o })) }))
agg <- merge(agg, adj, by = c("cell", "dep", "rel")); agg$R <- R
write.csv(agg, file.path(OUT, "E4d_rivals_summary.csv"), row.names = FALSE)
options(width = 250); cat(sprintf("\n=== E4d | rivals on identical data | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c(vars, paste0(vars, "_adj"))) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[order(sh$cell, sh$dep, sh$rel), c("cell", "dep", "rel", "grp", "mle_exists", "mle_hl", "mle_edge", "cvhl", "calslope", "spieg", "uss", "unc_hl", "unc_hl_adj", "calslope_adj", "spieg_adj", "uss_adj")], row.names = FALSE)
cat("\nrel = 0 rows are SIZE.  *_adj = size-adjusted power (critical value = empirical 5% null quantile of the same test in the same cell).\n")
