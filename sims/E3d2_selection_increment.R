# E3d2 -- the SELECTION increment as a within-model difference.
# Mean-field model (Bellec Thm 4.3): U,Z ~ N(0,1); y ~ Bern(expit(alpha + s U)); eta_hat = prox[gamma l_y](a U + sigma Z).
# "Random grouping"  = deciles of eta_hat (depends on y and Z).
# "Fixed grouping"   = deciles of the NOISE-FREE fitted index eta_bar(U) = prox[gamma l_{p1(U)}](a U), a function of U only.
# Model selection increment  Delta_g = E_MF[r_g | random] - E_MF[r_g | fixed].
# Observed target: the same difference from R replications (random on pi_hat(y) vs fixed on deciles of eta(beta_hat(pi_0))).
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 300; W <- if (length(args) > 1) as.integer(args[2]) else 6
CELL <- if (length(args) > 2) args[3] else "B416"
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
source(file.path(ROOT, "R", "closedform_engine.R"))
cfg <- switch(CELL, B416 = list(n = 400, p = 100, lambda = 416), B1000 = list(n = 400, p = 100, lambda = 1000), B50 = list(n = 400, p = 100, lambda = 50),
              P40 = list(n = 400, p = 40, lambda = 416), P160 = list(n = 400, p = 160, lambda = 416), A100 = list(n = 500, p = 5, lambda = 100))
n <- cfg$n; p <- cfg$p; lambda <- cfg$lambda; SEED <- 20260905
if (p == 5) { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
grp <- function(v, G = 10) pmin(ceiling(rank(v, ties.method = "first") / (length(v) / G)), G)
one <- function(r) {
  set.seed(SEED + r); X <- matrix(rnorm(n * p), n) %*% Sc; eta0 <- drop(X %*% b0); pi0 <- plogis(eta0); y <- rbinom(n, 1, pi0)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, p)))
  bbar <- cf_ridge(X1, pi0, lambda, Dpen); gfix <- grp(drop(X1 %*% bbar))
  pr <- cf_pieces(X, y, lambda); pf <- cf_pieces(X, y, lambda, groups = gfix)
  b <- cf_bellec_adjustments(X, y, pr); dn <- cf_denoised_w0(X, y, pr)
  c(a = sqrt(b$a2), sigma = sqrt(b$s2), gamma = b$gamma_hat, s = dn$c * sqrt(b$a2), alpha = dn$alpha, mean_eta = b$mean_eta,
    r_rand = pr$r, r_fix = pf$r, v_rand = pr$v, v_fix = pf$v)
}
cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, parLapply(cl, seq_len(R), one)); stopCluster(cl)
med <- apply(res, 2, median); a <- med["a"]; sigma <- med["sigma"]; gamma <- med["gamma"]; s <- med["s"]; alpha <- med["alpha"]; me <- med["mean_eta"]
# ---- mean-field quadrature, both groupings on the same nodes ----------------------------------
gh <- statmod_gh(64)
prox_logit <- function(t, y, g) { u <- t; for (it in 1:60) { e <- plogis(u); u <- u - (u - t + g * (e - y)) / (1 + g * e * (1 - e)) }; u }
grid <- expand.grid(iu = seq_along(gh$x), iz = seq_along(gh$x), y = c(0, 1))
U <- gh$x[grid$iu]; Z <- gh$x[grid$iz]; yv <- grid$y
p1 <- plogis(alpha + s * U); wgt <- gh$w[grid$iu] * gh$w[grid$iz] * ifelse(yv == 1, p1, 1 - p1)
eta_hat <- prox_logit(a * U + sigma * Z, yv, gamma) + me
eta_bar <- prox_logit(a * U, p1, gamma) + me                 # noise-free fitted index: function of U only
pi_hat <- plogis(eta_hat); resid <- yv - pi_hat; wvar <- pi_hat * (1 - pi_hat)
groups_of <- function(key) { o <- order(key); cw <- cumsum(wgt[o]) / sum(wgt); g <- integer(length(key)); g[o] <- pmin(findInterval(cw, seq(.1, .9, .1)) + 1, 10); g }
Er <- function(g) { Pg <- tapply(wgt, g, sum) / sum(wgt); mg <- tapply(wgt * resid, g, sum) / tapply(wgt, g, sum); vg <- tapply(wgt * wvar, g, sum) / tapply(wgt, g, sum); as.numeric(sqrt(n * Pg) * mg / sqrt(vg)) }
Er_rand <- Er(groups_of(eta_hat)); Er_fix <- Er(groups_of(eta_bar))
obs <- colMeans(res); o_rr <- obs[grep("^r_rand", names(obs))]; o_rf <- obs[grep("^r_fix", names(obs))]; o_vr <- obs[grep("^v_rand", names(obs))]; o_vf <- obs[grep("^v_fix", names(obs))]
D_model <- Er_rand - Er_fix; D_obs_r <- o_rr - o_rf; D_obs_v <- o_vr - o_vf
cat(sprintf("\n=== E3d2 | %s (n=%d p=%d lam=%g) | R = %d | medians: a=%.3f sigma=%.3f gamma=%.3f s=%.3f ===\n", CELL, n, p, lambda, R, a, sigma, gamma, s))
cat("model  E[r_g] random :", sprintf("%+.2f", Er_rand), "\nmodel  E[r_g] fixed  :", sprintf("%+.2f", Er_fix), "\n")
cat("obs    mean r random :", sprintf("%+.2f", o_rr), "\nobs    mean r fixed  :", sprintf("%+.2f", o_rf), "\n")
cat("\nSELECTION INCREMENT  model (rand-fix):", sprintf("%+.3f", D_model), "\n                     obs r (rand-fix):", sprintf("%+.3f", D_obs_r), "\n                     obs v (rand-fix):", sprintf("%+.3f", D_obs_v), "\n")
cat(sprintf("\n||model D||^2 = %.3f   ||obs D_r||^2 = %.3f   ||obs D_v||^2 = %.3f   corr(model D, obs D_v) = %.3f   corr(model D, obs D_r) = %.3f\n",
            sum(D_model^2), sum(D_obs_r^2), sum(D_obs_v^2), cor(D_model, D_obs_v), cor(D_model, D_obs_r)))
cat(sprintf("||obs mean v random||^2 = %.3f   ||obs mean v fixed||^2 = %.3f   ||obs mean v random - model D||^2 = %.3f  (what a corrected reference would leave)\n",
            sum(o_vr^2), sum(o_vf^2), sum((o_vr - D_model)^2)))
write.csv(data.frame(g = 1:10, Er_rand, Er_fix, D_model, o_rr, o_rf, o_vr, o_vf), file.path(ROOT, "results", sprintf("E3d2_selection_%s.csv", CELL)), row.names = FALSE)
