# E3d -- Does the Bellec/Sur-Candes MEAN-FIELD model predict the per-group residual means?
# Model (Thm 4.3, logistic): U,Z ~ N(0,1) indep; y ~ Bern(expit(alpha + s U)); eta_hat = prox[gamma l_y](a U + sigma Z).
# Inputs (a, sigma, gamma, s, alpha) are OBSERVABLE (Bellec 3.20 + the 1-D scale fit).  Prediction by 2-D
# Gauss-Hermite quadrature of E[r_g] = sqrt(n P_g) E[y - expit(eta_hat) | g] / sqrt(E[w | g]) for deciles g of eta_hat.
# Compared with the OBSERVED mean of r_g (and of v_g = r_g - mu_hat_g) over R replications, design B lambda = 416.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 300; W <- if (length(args) > 1) as.integer(args[2]) else 6
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE)
source(file.path(ROOT, "R", "closedform_engine.R"))
p <- 100; n <- 400; lambda <- 416; SEED <- 20260905
Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2))
one <- function(r) {
  set.seed(SEED + r); X <- matrix(rnorm(n * p), n) %*% Sc; eta0 <- drop(X %*% b0); y <- rbinom(n, 1, plogis(eta0))
  pc <- cf_pieces(X, y, lambda); b <- cf_bellec_adjustments(X, y, pc); dn <- cf_denoised_w0(X, y, pc)
  c(a = sqrt(b$a2), sigma = sqrt(b$s2), gamma = b$gamma_hat, s = dn$c * sqrt(b$a2), alpha = dn$alpha, mean_eta = b$mean_eta,
    true_s = sd(eta0), r = pc$r, v = pc$v, q = quantile(drop(cbind(1, X) %*% pc$beta), probs = seq(0.1, 0.9, 0.1)))
}
cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, parLapply(cl, seq_len(R), one)); stopCluster(cl)
obs <- colMeans(res); med <- apply(res, 2, median)
a <- med["a"]; sigma <- med["sigma"]; gamma <- med["gamma"]; s <- med["s"]; alpha <- med["alpha"]   # medians: robust to any divergent scale fit
cat(sprintf("scale estimate s: median %.3f, mean %.3f, 10%%-90%% quantiles %.3f-%.3f
", med["s"], obs["s"], quantile(res[, "s"], .1), quantile(res[, "s"], .9)))
cat(sprintf("\n=== E3d | B lam416 | R = %d ===\nobservables (mean over reps): a = %.3f  sigma = %.3f  gamma = %.3f  s = %.3f (true sd eta0 %.3f)  alpha = %.3f\n",
            R, a, sigma, gamma, s, obs["true_s"], alpha))

# ---- mean-field prediction by quadrature -------------------------------------------------
gh <- statmod_gh(48)
prox_logit <- function(t, y, g) { u <- t; for (it in 1:50) { e <- plogis(u); f <- u - t + g * (e - y); u <- u - f / (1 + g * e * (1 - e)) }; u }
grid <- expand.grid(iu = seq_along(gh$x), iz = seq_along(gh$x), y = c(0, 1))
U <- gh$x[grid$iu]; Z <- gh$x[grid$iz]; yv <- grid$y
p1 <- plogis(alpha + s * U); wgt <- gh$w[grid$iu] * gh$w[grid$iz] * ifelse(yv == 1, p1, 1 - p1)
eta_hat <- prox_logit(a * U + sigma * Z, yv, gamma) + obs["mean_eta"]
pi_hat <- plogis(eta_hat); resid <- yv - pi_hat; wvar <- pi_hat * (1 - pi_hat)
o <- order(eta_hat); cw <- cumsum(wgt[o]) / sum(wgt); grp <- pmin(findInterval(cw, seq(0.1, 0.9, 0.1)) + 1, 10)
g_of <- integer(length(eta_hat)); g_of[o] <- grp
Pg <- tapply(wgt, g_of, sum) / sum(wgt); mg <- tapply(wgt * resid, g_of, sum) / tapply(wgt, g_of, sum); vg <- tapply(wgt * wvar, g_of, sum) / tapply(wgt, g_of, sum)
pred_r <- sqrt(n * Pg) * mg / sqrt(vg)
qpred <- sapply(seq(0.1, 0.9, 0.1), function(q) eta_hat[o][which(cw >= q)[1]])
obs_r <- obs[grep("^r", names(obs))]; obs_v <- obs[grep("^v", names(obs))]; obs_q <- obs[grep("^q", names(obs))]
cat("\nmodel deciles of eta_hat :", sprintf("%+.2f", qpred), "\nobserved deciles         :", sprintf("%+.2f", obs_q), "\n")
cat("\npredicted  E[r_g] (mean-field):", sprintf("%+.3f", pred_r), "\nobserved   mean r_g            :", sprintf("%+.3f", obs_r),
    "\nobserved   mean v_g = r - mu_hat:", sprintf("%+.3f", obs_v), "\n")
cat(sprintf("\n||pred E r||^2 = %.3f   ||obs mean r||^2 = %.3f   ||obs mean v||^2 = %.3f   corr(pred, obs r) = %.3f\n",
            sum(pred_r^2), sum(obs_r^2), sum(obs_v^2), cor(pred_r, obs_r)))
write.csv(data.frame(group = 1:10, pred_r = as.numeric(pred_r), obs_r = as.numeric(obs_r), obs_v = as.numeric(obs_v)), file.path(ROOT, "results", "E3d_meanfield_group_means.csv"), row.names = FALSE)
