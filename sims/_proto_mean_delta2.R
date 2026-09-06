# Prototype: second-order delta-method estimate of the non-centrality delta = E[v] under the null world,
# deterministic (no random numbers).  y* = pi0 + W0^{1/2} eps; beta_hat(y*) ~= beta_bar + M^{-1} X'(y* - pi0);
# E[v] ~= v(pi0) + (1/2) sum_i w0_i d^2 v / d y_i^2   (central finite differences, 2n cheap evaluations of f).
source("R/closedform_engine.R")
set.seed(3); p <- 100; n <- 400
Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2))
X <- matrix(rnorm(n * p), n) %*% Sc; eta0 <- drop(X %*% b0); pi0 <- plogis(eta0); y <- rbinom(n, 1, pi0)
lambda <- 416
pc <- cf_pieces(X, y, lambda); X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, p))); idx <- split(seq_len(n), pc$groups)

mean_v_delta2 <- function(pi_gen, w_gen, h = 0.25) {
  bbar <- cf_ridge(X1, pi_gen, lambda, Dpen)                       # fit on the noise-free response
  Mi   <- solve(crossprod(X1, pmax(plogis(drop(X1 %*% bbar)) * (1 - plogis(drop(X1 %*% bbar))), 1e-8) * X1) + lambda * Dpen)
  f <- function(ystar) cf_v_of_beta(ystar, X1, bbar + drop(Mi %*% crossprod(X1, ystar - pi_gen)), lambda, Dpen, idx)
  v0 <- f(pi_gen); acc <- 0
  for (i in seq_len(n)) { e <- numeric(n); e[i] <- h; acc <- acc + w_gen[i] * (f(pi_gen + e) - 2 * v0 + f(pi_gen - e)) / h^2 }
  list(v0 = v0, delta = v0 + 0.5 * acc)
}
# generator 1: the TRUE pi0 (oracle, to validate the method); generator 2: pi(beta_tilde); generator 3: de-noised pi0-hat
pi_t <- pmin(pmax(plogis(drop(X1 %*% pc$beta_tilde)), 1e-6), 1 - 1e-6)
dn <- cf_denoised_w0(X, y, pc); b <- cf_bellec_adjustments(X, y, pc); k <- b$a2 / (b$a2 + b$s2); muz <- b$m * k; tau <- sqrt(b$a2 * b$s2 / (b$a2 + b$s2))
gh <- statmod_gh(20); pi5 <- numeric(n); for (q in 1:20) pi5 <- pi5 + gh$w[q] * plogis(dn$alpha + dn$c * (muz + tau * gh$x[q]))
for (nm in c("true pi0", "pi(beta_tilde)", "de-noised pi0-hat")) {
  pg <- switch(nm, "true pi0" = pi0, "pi(beta_tilde)" = pi_t, "de-noised pi0-hat" = pi5)
  t0 <- Sys.time(); r <- mean_v_delta2(pg, pg * (1 - pg))
  cat(sprintf("%-18s ||v(y*)||^2 = %.3f   ||delta2||^2 = %.3f   (E3 measured ||E v||^2 ~ 0.28 at this cell)   %.1fs\n", nm, sum(r$v0^2), sum(r$delta^2), as.numeric(Sys.time() - t0)))
}
