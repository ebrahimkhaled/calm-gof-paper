# =====================================================================
# closedform_engine.R -- the R1 engine
#
# Exact (no-bootstrap) reference for the shrinkage-corrected grouped
# statistics.  Under H0^str, to first order,
#
#     r - mu_hat  ~  N(0, Omega_MLE),      Omega_MLE = I_G - U F^{-1} U'
#
# (Series C paper, Proposition 2).  Hence
#     SC.HL   = ||r - mu_hat||^2                  ~ sum_k  lambda_k chi^2_1,   lambda_k = eig(Omega_MLE)
#     SC.EDGE = (r - mu_hat)' P_Z (r - mu_hat)    ~ sum_k  nu_k     chi^2_1,   nu_k     = eig(Omega^{1/2} P_Z Omega^{1/2})
# and the p-value is the exact tail of a weighted chi-square (Davies 1980,
# via CompQuadForm), with NO chi^2_{G-2} approximation anywhere.
#
# `deflate` multiplies Omega_MLE.  deflate = 1 is the first-order reference
# (row (a) of Table S1 in the Series C paper: valid at fixed p, conservative
# at p/n = 0.25).  The proportional-regime adjustment will enter here.
#
# The engine also returns the OBSERVABLE quantities of Bellec (2022,
# arXiv 2204.06990, eq. 3.20) in our parameterization, for the deflation
# diagnostic:  df_M = tr[X M^{-1} X' W],  trV_M = tr(W) - tr[W X M^{-1} X' W],
# v_hat = trV_M / n,  r2_hat = ||psi||^2 / n,  gamma_hat = df_M / trV_M,
# and the same with F in place of M (the unpenalized analogues).
#
# lambda is on the THEORY scale: lambda = n * lambda_glmnet.
# =====================================================================
suppressPackageStartupMessages(library(CompQuadForm))

cf_expit <- function(z) 1 / (1 + exp(-z))

cf_ridge <- function(X1, y, lambda, Dpen, tol = 1e-10, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    z  <- eta + (y - pi) / w
    bn <- drop(solve(crossprod(X1, w * X1) + lambda * Dpen, crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  beta
}

cf_qf <- function(u, Z) { zr <- crossprod(Z, u); as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr) }

# ---- everything the tests need, from one fit ----------------------------
# groups: optional integer vector of group labels (to group on the TRUE index,
#         or on the debiased index, in diagnostics).  Default: deciles of pi_hat.
cf_pieces <- function(X, y, lambda, G = 10, beta = NULL, groups = NULL) {
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, ncol(X))))
  if (is.null(beta)) beta <- cf_ridge(X1, y, lambda, Dpen)
  eta <- drop(X1 %*% beta)
  pi  <- pmin(pmax(cf_expit(eta), 1e-6), 1 - 1e-6)
  w   <- pmax(pi * (1 - pi), 1e-8)
  if (is.null(groups)) groups <- pmin(ceiling(rank(pi, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg  <- vapply(idx, function(I) sum(w[I]), 0.0)
  r   <- (vapply(idx, function(I) sum(y[I]), 0.0) - vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U   <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]), numeric(ncol(X1)))) / sqrt(Vg)
  Fm  <- crossprod(X1, w * X1); K <- lambda * Dpen; M <- Fm + K
  Fi  <- solve(Fm); Mi <- solve(M)
  bt  <- beta + drop(Fi %*% (K %*% beta))                 # one-step debias
  mu  <- drop(U %*% Mi %*% (K %*% bt))                    # the shrinkage displacement
  v   <- r - mu                                           # the corrected residual
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0)
  Z   <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  Om_MLE <- diag(G) - U %*% Fi %*% t(U)
  Om_K   <- diag(G) - U %*% Mi %*% (Fm + 2 * K) %*% Mi %*% t(U)

  # ---- Bellec-type observables (traces computed without any n x n matrix)
  XW2X  <- crossprod(X1, w^2 * X1)
  df_M  <- sum(diag(Mi %*% Fm));  trV_M <- sum(w) - sum(diag(Mi %*% XW2X))
  df_F  <- ncol(X1);              trV_F <- sum(w) - sum(diag(Fi %*% XW2X))
  psi   <- y - pi
  obs <- c(n = n, p = ncol(X), G = G, lambda = lambda,
           df_M = df_M, trV_M = trV_M, v_hat = trV_M / n, r2_hat = sum(psi^2) / n,
           gamma_hat = df_M / trV_M, df_F = df_F, trV_F = trV_F,
           trW = sum(w), tr_OmMLE = sum(diag(Om_MLE)), tr_OmK = sum(diag(Om_K)),
           norm_bhat = sqrt(sum(beta[-1]^2)), norm_btilde = sqrt(sum(bt[-1]^2)))

  list(S_dec = sum(v^2), S_edge = if (is.null(Z)) NA_real_ else cf_qf(v, Z),
       S_dec_unc = sum(r^2), S_edge_unc = if (is.null(Z)) NA_real_ else cf_qf(r, Z),
       r = r, mu = mu, v = v, Om_MLE = Om_MLE, Om_K = Om_K, Z = Z,
       beta = beta, beta_tilde = bt, pi = pi, groups = groups, obs = obs)
}

# ---- exact tail probability of a Gaussian quadratic form ----------------
# S = u' A u with u ~ N(0, Sigma).  A = I (decile) or P_Z (EDGE).
cf_pvalue <- function(S, Sigma, Z = NULL) {
  if (is.na(S)) return(NA_real_)
  if (is.null(Z)) {
    lam <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  } else {
    e  <- eigen(Sigma, symmetric = TRUE)
    Sh <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))      # Sigma^{1/2}
    P  <- Z %*% solve(crossprod(Z), t(Z))
    lam <- eigen(Sh %*% P %*% Sh, symmetric = TRUE, only.values = TRUE)$values
  }
  lam <- lam[lam > 1e-9]
  if (!length(lam)) return(NA_real_)
  d <- CompQuadForm::davies(S, lambda = lam, lim = 1e5, acc = 1e-7)
  p <- d$Qq
  if (d$ifault != 0 || is.na(p) || p < 0 || p > 1) p <- CompQuadForm::imhof(S, lambda = lam)$Qq
  min(max(p, 0), 1)
}

# ---- the closed-form test ----------------------------------------------
# FINAL DESIGN (E9/E9b, 2026-09-07): reference R5 + basis-specific kappa-inflation; bases = decile (SC.HL) and the
# attenuation-adapted EDGE (SC.EDGE.adaptive).  Valid (size .016-.060) in all eight registered cells, kappa up to 0.40.
# CALM -- Calibration Assessment under Lambda-shrunk Models (named 2026-09-07).  The exact analytic reference for the
# shrinkage-corrected statistics SC.HL and SC.EDGE.  closedform.gof() is kept as an alias.
calm.gof <- function(X, y, lambda, G = 10, basis = c("decile", "adaptive", "edge"),
                     lambda_scale = c("theory", "glmnet"),   # "glmnet": lambda is glmnet's (1/n scale) -> converted to n * lambda; BOTH printed
                           reference = c("R5", "R2", "R4", "R0", "R1", "R3"),
                           grouping = c("fitted", "heavy", "loo"),   # "heavy" = deciles of the index of a 5x-heavier ridge fit (E5: removes the selection term, no constants); "loo" is INVALID (E3e), kept for reproducibility
                           heavy_mult = 5,
                           inflate = c("kappa", "none", "hbar", "hr"),
                           degree_rule = c("rho2k", "rho2", "kappa", "rho"),   # EDGE degree: "rho2" (default) = keep the cubic iff rho_hat^6 >= .2, floor at 2; "kappa" = 3 at fixed p else 2; "rho" = retired rule (can drop to 1)
                           deflate = 1, groups = NULL, uncorrected = FALSE) {
  basis <- match.arg(basis, c("decile", "adaptive", "edge"), several.ok = TRUE)
  reference <- match.arg(reference); grouping <- match.arg(grouping); lambda_scale <- match.arg(lambda_scale)
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  lambda_glmnet <- if (lambda_scale == "glmnet") lambda else lambda / n
  if (lambda_scale == "glmnet") lambda <- n * lambda          # the #1 trap in this project: glmnet's lambda is on the 1/n scale
  if (is.null(groups) && grouping == "heavy") {
    X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, ncol(X))))
    groups <- pmin(ceiling(rank(drop(X1 %*% cf_ridge(X1, y, heavy_mult * lambda, Dpen)), ties.method = "first") / (n / G)), G)
  }
  if (is.null(groups) && grouping == "loo") {
    # one-step leave-one-out fitted index (Pregibon 1981): the group of i does not depend on y_i
    X1 <- cbind(1, X); Dpen <- diag(c(0, rep(1, ncol(X))))
    beta <- cf_ridge(X1, y, lambda, Dpen); eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    h <- w * rowSums((X1 %*% solve(crossprod(X1, w * X1) + lambda * Dpen)) * X1)
    eta_loo <- eta - h * (y - pi) / (w * pmax(1 - h, 1e-6))
    groups <- pmin(ceiling(rank(eta_loo, ties.method = "first") / (n / G)), G)
  }
  pc  <- cf_pieces(X, y, lambda, G = G, groups = groups)
  # reference covariance for v = r - mu_hat:  R0 = Omega_MLE (first order); R1-R4 cf_references(); R5 de-noised (Bellec)
  Sig <- if (reference == "R0") pc$Om_MLE else if (reference == "R5") cf_reference_R5(X, y, lambda, pc)$R5 else cf_references(X, y, lambda, pc)[[reference]]
  # kappa-inflation (2026-09-07): the selection non-centrality is O(kappa); inflating the reference by (1 + g_hat)
  # matches E[S] and errs toward conservatism (scaled central law has the heavier tail). Laws fitted in sims/_delta_law.R.
  # The selection non-centrality is an S-shape that is almost entirely CUBIC in the group index (E9/E9b, 2026-09-07):
  # its fraction of the reference trace is ~0.24 kappa for the decile basis, ~1.15 kappa for EDGE-3, ~0.11 kappa for
  # EDGE degree <= 2.  Inflation is therefore BASIS-SPECIFIC.  (hbar / hr variants kept for the record.)
  inflate <- match.arg(inflate); kappa <- ncol(X) / n; hbar <- unname(pc$obs["df_M"]) / n
  g_of <- function(deg) if (inflate == "none") 0 else switch(inflate,
    kappa = kappa * c(dec = 0.24, e1 = 0.11, e2 = 0.11, e3 = 1.15)[[deg]],
    hbar  = 0.28 * sqrt(hbar) * c(dec = 1, e1 = 0.5, e2 = 0.5, e3 = 4.8)[[deg]],
    hr    = { b <- cf_bellec_adjustments(X, y, pc); (0.50 * hbar + 0.07 * (1 - b$a2 / (b$a2 + b$s2))) * c(dec = 1, e1 = 0.5, e2 = 0.5, e3 = 4.8)[[deg]] })
  Sig0 <- deflate * Sig; Sig <- (1 + g_of("dec")) * Sig0; g_hat <- g_of("dec")
  out <- list(lambda = lambda, lambda_glmnet = lambda_glmnet, G = G, reference = reference, grouping = grouping, inflate = inflate, g_hat = g_hat, deflate = deflate, obs = pc$obs)
  if ("decile" %in% basis) {
    out$SC.HL <- list(statistic = pc$S_dec, p.value = cf_pvalue(pc$S_dec, Sig))
    if (uncorrected) out$SC.HL$p.uncorrected <- cf_pvalue(pc$S_dec_unc, pc$Om_MLE)
  }
  if ("edge" %in% basis) {
    out$SC.EDGE <- list(statistic = pc$S_edge, p.value = cf_pvalue(pc$S_edge, (1 + g_of("e3")) * Sig0, pc$Z), degree = 3L, g_hat = g_of("e3"))
    if (uncorrected) out$SC.EDGE$p.uncorrected <- cf_pvalue(pc$S_edge_unc, pc$Om_MLE, pc$Z)
  }
  if ("adaptive" %in% basis) {
    # Degree rule (2026-09-04, E4c/E4c2): the rho_hat rule is DEAD as a default -- under a misspecified model rho_hat collapses
    # (Bellec's a^2, sigma^2 are null quantities), the rule dropped to degree 1 and lost 25-50 pp of power against a quadratic
    # departure. The degree must depend on the DESIGN only: k = 3 at fixed p (p/n < 0.05), k = 2 in the proportional regime,
    # where the cubic direction carries ~90% of the selection non-centrality and rho^6 ~ 0.1 of any cubic signal (E9b, E8).
    # E4c2 on identical data: EDGE-2 holds size (.018-.042) and matches/beats prepivot at every kappa. rho_hat is still REPORTED.
    # FINAL RULE (E4c2 + E4c3 + glaucoma, 2026-09-04 17:00): "rho2" = the rho_hat rule with a FLOOR at degree 2 -- keep the cubic only
    # when rho_hat^6 >= 0.2 (rho_hat >= 0.765), never drop to degree 1. The kappa-only rule kept the level and the power but LOST the
    # glaucoma misfit (EDGE-2 p = .57 vs EDGE-3 p = .033 at rho_hat = .92): a design-only rule ignores signal strength. Offline check on
    # identical data: size .018-.050 in all cells; power = EDGE-2 at kappa (.42/.64 vs prepivot .46/.62 at kappa .25), EDGE-3 at fixed p.
    # FINAL RULE (2026-09-05, after E8 exposed the remaining hole): keep the cubic when EITHER there is no selection term to avoid
    # (kappa < 0.05: the term is O(kappa) and negligible) OR the index is accurate enough to carry a cubic signal (rho_hat^6 >= tau).
    # Floor at 2 always. The kappa clause matters because rho_hat is a NULL quantity: under a strong cubic departure at fixed p it
    # falls from .94 to .75 (E8), crossing the threshold and discarding the very direction that carries the signal (power .68 -> .08).
    # The rho clause matters because a design-only rule drops the cubic on the glaucoma data, where rho_hat = .92 and the misfit IS cubic.
    degree_rule <- match.arg(degree_rule)
    ad <- cf_adaptive_degree(X, y, pc)
    kdeg <- switch(degree_rule,
                   rho2k = if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L,
                   rho2  = max(2L, ad$degree),
                   kappa = if (kappa < 0.05) 3L else 2L,
                   rho   = ad$degree)
    Zk <- cf_edge_basis(pc, kdeg); gk <- g_of(c("e1", "e2", "e3")[kdeg])
    Sk <- if (is.null(Zk)) NA_real_ else cf_qf(pc$v, Zk)
    out$SC.EDGE.adaptive <- list(statistic = Sk, p.value = cf_pvalue(Sk, (1 + gk) * Sig0, Zk), degree = kdeg, rho_hat = ad$rho_hat, g_hat = gk, degree_rule = degree_rule)
  }
  class(out) <- "calm.gof"; out
}
closedform.gof <- calm.gof   # alias

print.calm.gof <- function(x, ...) {
  cat("\nCALM: exact reference for the shrinkage-corrected Hosmer-Lemeshow test (SC.HL, SC.EDGE)\n")
  cat(sprintf("n = %d, p = %d, lambda = %.4g (theory scale) = %.4g (glmnet scale), G = %d
reference = %s, grouping = %s, inflate = %s (g = %.3f)

", x$obs["n"], x$obs["p"], x$lambda, x$lambda_glmnet, x$G, x$reference, x$grouping, x$inflate, x$g_hat))
  for (nm in c("SC.HL", "SC.EDGE", "SC.EDGE.adaptive")) if (!is.null(x[[nm]])) {
    cat(sprintf("  %-16s statistic = %8.4f   p = %.4f", nm, x[[nm]]$statistic, x[[nm]]$p.value))
    if (!is.null(x[[nm]]$degree) && nm == "SC.EDGE.adaptive") cat(sprintf("   [degree %d, rho_hat %.2f]", x[[nm]]$degree, x[[nm]]$rho_hat))
    if (!is.null(x[[nm]]$p.uncorrected)) cat(sprintf("   (uncorrected p = %.5f)", x[[nm]]$p.uncorrected))
    cat("\n")
  }
  cat("\n"); invisible(x)
}

# =====================================================================
# E3 additions -- OBSERVABLE, matrix-level references for v = r - mu_hat
#
# Two ingredients, used factorially:
#   (J)  the EXACT Jacobian of v with respect to y at the observed data, by
#        implicit differentiation:  d beta_hat / d y = M^{-1} X'  (from the
#        penalized score), and d f / d beta by central finite differences over
#        the p+1 coordinates, where f(y, beta) = r - mu_hat with the GROUPING
#        HELD FIXED.  The first-order A* = V^{-1/2} C - U F^{-1} X' drops the
#        dependence of W, V_g, F and mu_hat on beta; J keeps all of it.
#   (W~) the Bernoulli variance evaluated at the DEBIASED generator pi(beta_tilde)
#        instead of the shrunken fit pi(beta_hat) -- the analytic analogue of
#        prepivoting's generator choice (Series C, Table S1 row (d)).
# Reference covariance = J W J'  with J in {A*, J_exact}, W in {W(beta_hat), W(beta_tilde)}.
# =====================================================================

# v = r - mu_hat as a function of (y, beta), grouping fixed
cf_v_of_beta <- function(y, X1, beta, lambda, Dpen, idx) {
  eta <- drop(X1 %*% beta); pi <- pmin(pmax(cf_expit(eta), 1e-6), 1 - 1e-6); w <- pmax(pi * (1 - pi), 1e-8)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r  <- (vapply(idx, function(I) sum(y[I]), 0.0) - vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U  <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]), numeric(ncol(X1)))) / sqrt(Vg)
  Fm <- crossprod(X1, w * X1); K <- lambda * Dpen; M <- Fm + K
  bt <- beta + drop(solve(Fm, K %*% beta))
  mu <- drop(U %*% solve(M, K %*% bt))
  r - mu
}

# all four reference covariances from one fit (pc = cf_pieces output)
cf_references <- function(X, y, lambda, pc, h = 1e-4) {
  X1 <- cbind(1, as.matrix(X)); n <- nrow(X1); q <- ncol(X1)
  Dpen <- diag(c(0, rep(1, q - 1))); idx <- split(seq_len(n), pc$groups); G <- length(idx)
  beta <- pc$beta
  pi_h <- pc$pi; w_h <- pmax(pi_h * (1 - pi_h), 1e-8)
  pi_t <- pmin(pmax(cf_expit(drop(X1 %*% pc$beta_tilde)), 1e-6), 1 - 1e-6); w_t <- pmax(pi_t * (1 - pi_t), 1e-8)
  Vg <- vapply(idx, function(I) sum(w_h[I]), 0.0)
  Cmat <- matrix(0, G, n); for (g in seq_len(G)) Cmat[g, idx[[g]]] <- 1
  VC <- Cmat / sqrt(Vg)                                            # V^{-1/2} C   (G x n)
  U  <- t(vapply(idx, function(I) colSums(w_h[I] * X1[I, , drop = FALSE]), numeric(q))) / sqrt(Vg)
  Fm <- crossprod(X1, w_h * X1); M <- Fm + lambda * Dpen
  Astar <- VC - U %*% solve(Fm, t(X1))                            # first-order Jacobian (G x n)
  # exact Jacobian: d f/d y = VC (beta fixed);  d f/d beta by central FD;  d beta/d y = M^{-1} X'
  dfdb <- matrix(0, G, q)
  for (j in seq_len(q)) {
    e <- rep(0, q); e[j] <- h
    dfdb[, j] <- (cf_v_of_beta(y, X1, beta + e, lambda, Dpen, idx) - cf_v_of_beta(y, X1, beta - e, lambda, Dpen, idx)) / (2 * h)
  }
  Jex <- VC + dfdb %*% solve(M, t(X1))                             # (G x n)
  # (R4) variance at the RECALIBRATED scale: fit y ~ a + b * eta_hat in-sample (2 parameters),
  #      then w_r = pi_r (1 - pi_r) with pi_r = expit(a + b eta_hat).  Repairs the shrinkage
  #      scale (b ~ 1/alpha) that inflates w_h toward 1/4, without inheriting beta_tilde's noise.
  eta_h <- drop(X1 %*% beta)
  rc <- tryCatch(coef(glm(y ~ eta_h, family = binomial)), error = function(e) c(0, 1))
  pi_r <- pmin(pmax(cf_expit(rc[1] + rc[2] * eta_h), 1e-6), 1 - 1e-6); w_r <- pmax(pi_r * (1 - pi_r), 1e-8)
  list(R0 = Astar %*% (w_h * t(Astar)),     # = Omega_MLE (identity check)
       R1 = Jex   %*% (w_h * t(Jex)),       # exact Jacobian, fitted variance
       R2 = Astar %*% (w_t * t(Astar)),     # first-order Jacobian, debiased-generator variance
       R3 = Jex   %*% (w_t * t(Jex)),       # both
       R4 = Astar %*% (w_r * t(Astar)),     # first-order Jacobian, RECALIBRATED-scale variance
       w_ratio = sum(w_t) / sum(w_h),       # generator variance / fitted variance
       w_ratio_recal = sum(w_r) / sum(w_h), # recalibrated variance / fitted variance
       recal_slope = unname(rc[2]))
}

# implied first-order mean of S under N(0, Sigma): tr(Sigma) for decile, tr(Sigma^{1/2} P_Z Sigma^{1/2}) for EDGE
cf_implied_df <- function(Sigma, Z = NULL) {
  if (is.null(Z)) return(sum(diag(Sigma)))
  e <- eigen(Sigma, symmetric = TRUE); Sh <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  P <- Z %*% solve(crossprod(Z), t(Z)); sum(diag(Sh %*% P %*% Sh))
}

# =====================================================================
# R5 -- the DE-NOISED null variance, from Bellec (2022) Theorem 4.3
#
# Thm 4.3 (logistic):  m_i = eta_hat_i - gamma* psi_hat_i  ~=  a* U_i + sigma* Z_i,
# U_i = the true (unit-variance) index, Z_i ~ N(0,1) INDEPENDENT.  With Gaussian
# design, (a* U | m) ~ N( m a^2/(a^2+s^2),  a^2 s^2/(a^2+s^2) ).  The true linear
# predictor is eta_0 = c * (a* U) for ONE unknown scale c (the signal strength),
# which the KNOWN logistic link identifies: fit c by 1-D maximum likelihood of y
# on the de-noised index (Gauss-Hermite).  Then
#     w0_i = E[ expit(c zeta)(1 - expit(c zeta)) | m_i ]
# is the observable estimate of the null Bernoulli variance pi_0(1-pi_0), and the
# reference is  R5 = A* W0 A*'.   Inputs: (gamma_hat, a_hat^2, sigma_hat^2) from
# Bellec eq. (3.20), computed with Sigma_hat = cov(X)  (p < n).
# =====================================================================
cf_bellec_adjustments <- function(X, y, pc, Sigma_x = NULL) {   # Sigma_x: pass the TRUE covariance for the oracle check (referee item)
  X  <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)                  # Bellec's model has no intercept
  eta <- drop(cbind(1, X) %*% pc$beta); psi <- y - pc$pi
  Xb  <- eta - mean(eta)                                        # X beta_hat, centred
  gamma_hat <- unname(pc$obs["gamma_hat"]); v_hat <- unname(pc$obs["v_hat"])
  r2 <- sum(psi^2) / n
  S  <- if (is.null(Sigma_x)) crossprod(Xc) / n else Sigma_x; eS <- eigen(S, symmetric = TRUE)
  Sinvhalf <- eS$vectors %*% (1 / sqrt(pmax(eS$values, 1e-8)) * t(eS$vectors))
  t2 <- sum((Sinvhalf %*% crossprod(Xc, psi))^2) / n^2 + (2 * v_hat / n) * sum(psi * Xb) +
        (v_hat^2 / n) * sum((Xb - gamma_hat * psi)^2) - (p / n) * r2
  a2 <- ((v_hat / n) * sum((Xb - gamma_hat * psi)^2) + sum(psi * Xb) / n - gamma_hat * r2)^2 / max(t2, 1e-8)
  s2 <- max(sum((Xb - gamma_hat * psi)^2) / n - a2, 1e-8)
  list(m = Xb - gamma_hat * psi, gamma_hat = gamma_hat, t2 = t2, a2 = a2, s2 = s2, mean_eta = mean(eta))
}

cf_denoised_w0 <- function(X, y, pc, nodes = 20) {
  b  <- cf_bellec_adjustments(X, y, pc)
  k  <- b$a2 / (b$a2 + b$s2)
  mu <- b$m * k + 0                      # conditional mean of a*U given m (centred index)
  tau <- sqrt(b$a2 * b$s2 / (b$a2 + b$s2))
  gh <- statmod_gh(nodes)                # Gauss-Hermite for N(0,1)
  # eta_0 = alpha + c * zeta ; fit (alpha, c) by ML of y on the de-noised index
  nll <- function(par) {
    alpha <- par[1]; c <- exp(par[2])
    ll <- 0
    for (q in seq_len(nodes)) {
      z <- mu + tau * gh$x[q]; pr <- pmin(pmax(cf_expit(alpha + c * z), 1e-9), 1 - 1e-9)
      ll <- ll + gh$w[q] * ifelse(y == 1, pr, 1 - pr)
    }
    -sum(log(pmax(ll, 1e-300)))
  }
  # bounded: log c in [log 0.05, log 20]; an unbounded BFGS diverged on a few draws (c -> 1e9, a step-function link)
  op <- optim(c(b$mean_eta, 0), nll, method = "L-BFGS-B", lower = c(-10, log(0.05)), upper = c(10, log(20)))
  alpha <- op$par[1]; c <- exp(op$par[2])
  w0 <- numeric(length(y))
  for (q in seq_len(nodes)) { z <- mu + tau * gh$x[q]; pr <- cf_expit(alpha + c * z); w0 <- w0 + gh$w[q] * pr * (1 - pr) }
  list(w0 = pmax(w0, 1e-8), c = c, alpha = alpha, est_signal_sd = c * sqrt(b$a2), a2 = b$a2, s2 = b$s2, gamma_hat = b$gamma_hat)
}

# Gauss-Hermite nodes/weights for a standard normal (probabilists'), no dependency
statmod_gh <- function(n) {
  i <- seq_len(n - 1); J <- matrix(0, n, n); J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2)
}

cf_reference_R5 <- function(X, y, lambda, pc) {
  X1 <- cbind(1, as.matrix(X)); n <- nrow(X1); q <- ncol(X1)
  idx <- split(seq_len(n), pc$groups); G <- length(idx)
  w_h <- pmax(pc$pi * (1 - pc$pi), 1e-8); Vg <- vapply(idx, function(I) sum(w_h[I]), 0.0)
  Cmat <- matrix(0, G, n); for (g in seq_len(G)) Cmat[g, idx[[g]]] <- 1
  U  <- t(vapply(idx, function(I) colSums(w_h[I] * X1[I, , drop = FALSE]), numeric(q))) / sqrt(Vg)
  Fm <- crossprod(X1, w_h * X1)
  Astar <- Cmat / sqrt(Vg) - U %*% solve(Fm, t(X1))
  dn <- cf_denoised_w0(X, y, pc)
  list(R5 = Astar %*% (dn$w0 * t(Astar)), w_ratio_R5 = sum(dn$w0) / sum(w_h),
       est_signal_sd = dn$est_signal_sd, a2 = dn$a2, s2 = dn$s2, c = dn$c)
}

# =====================================================================
# Package item F -- GENERAL SMOOTH PENALTY (Series C Remark 2 made numerical)
# Estimating equation:  X'(y - pi(beta)) + a(beta) = 0,   K(beta) = -d a / d beta  (p+1 x p+1)
#   ridge : a(beta) = -lambda D beta,            K = lambda D
#   gridge: a(beta) = -diag(lambda_j) beta,      K = diag(lambda_j)    (generalized ridge)
#   firth : a_j(beta) = (1/2) tr(F^{-1} dF/dbeta_j),  dF/dbeta_j = X' diag(w (1-2pi) x_j) X   (Jeffreys score)
# Displacement mu = U M^{-1} a(beta_hat) sign-adjusted (Remark 2: with a = -K beta this is U M^{-1} K beta),
# debiased beta_tilde = beta_hat + F^{-1} K beta_hat -> general form beta_hat - F^{-1} a(beta_hat), M = F + K.
# =====================================================================
cf_penalty <- function(type = c("ridge", "gridge", "firth"), lambda = NULL, lambda_vec = NULL) {
  type <- match.arg(type)
  list(type = type,
       a = function(beta, X1, w, pi) switch(type,
         ridge  = -lambda * c(0, beta[-1]),
         gridge = -c(0, lambda_vec) * beta,
         firth  = { Fi <- solve(crossprod(X1, w * X1)); q <- ncol(X1)
                    vapply(seq_len(q), function(j) 0.5 * sum(Fi * crossprod(X1, (w * (1 - 2 * pi) * X1[, j]) * X1)), 0.0) }),
       K = function(beta, X1, w, pi, h = 1e-4) switch(type,
         ridge  = lambda * diag(c(0, rep(1, ncol(X1) - 1))),
         gridge = diag(c(0, lambda_vec)),
         firth  = { q <- ncol(X1); Km <- matrix(0, q, q)
                    for (j in seq_len(q)) { e <- numeric(q); e[j] <- h
                      ap <- { eta <- drop(X1 %*% (beta + e)); p1 <- cf_expit(eta); w1 <- pmax(p1 * (1 - p1), 1e-8); Fi <- solve(crossprod(X1, w1 * X1)); vapply(seq_len(q), function(k) 0.5 * sum(Fi * crossprod(X1, (w1 * (1 - 2 * p1) * X1[, k]) * X1)), 0.0) }
                      am <- { eta <- drop(X1 %*% (beta - e)); p1 <- cf_expit(eta); w1 <- pmax(p1 * (1 - p1), 1e-8); Fi <- solve(crossprod(X1, w1 * X1)); vapply(seq_len(q), function(k) 0.5 * sum(Fi * crossprod(X1, (w1 * (1 - 2 * p1) * X1[, k]) * X1)), 0.0) }
                      Km[, j] <- -(ap - am) / (2 * h) }
                    (Km + t(Km)) / 2 }))
}
# fit the penalized estimating equation by Newton with the penalty's own curvature
cf_fit_penalty <- function(X1, y, pen, tol = 1e-9, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- cf_expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    s <- drop(crossprod(X1, y - pi)) + pen$a(beta, X1, w, pi)
    # Newton curvature: F + K for the cheap penalties; for Firth use F alone (the Jeffreys curvature is O(1/n) and its
    # finite-difference evaluation is O(p^2) full passes -- computing it every iteration made E7 run for hours; logistf
    # iterates with F as well). K is still evaluated ONCE, at convergence, in cf_pieces_penalty().
    Mh <- crossprod(X1, w * X1) + if (pen$type == "firth") 0 else pen$K(beta, X1, w, pi)
    step <- drop(solve(Mh, s)); beta <- beta + step
    if (max(abs(step)) < tol) break
  }
  beta
}
cf_pieces_penalty <- function(X, y, pen, G = 10, groups = NULL) {
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y); X1 <- cbind(1, X)
  beta <- cf_fit_penalty(X1, y, pen)
  eta <- drop(X1 %*% beta); pi <- pmin(pmax(cf_expit(eta), 1e-6), 1 - 1e-6); w <- pmax(pi * (1 - pi), 1e-8)
  if (is.null(groups)) groups <- pmin(ceiling(rank(pi, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg <- vapply(idx, function(I) sum(w[I]), 0.0)
  r  <- (vapply(idx, function(I) sum(y[I]), 0.0) - vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U  <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]), numeric(ncol(X1)))) / sqrt(Vg)
  Fm <- crossprod(X1, w * X1); a <- pen$a(beta, X1, w, pi); K <- pen$K(beta, X1, w, pi); M <- Fm + K
  Fi <- solve(Fm)
  bt <- beta - drop(Fi %*% a)                      # general debias: beta_tilde = beta_hat - F^{-1} a(beta_hat)
  eta_t <- drop(X1 %*% bt); pi_t <- pmin(pmax(cf_expit(eta_t), 1e-6), 1 - 1e-6)
  a_t <- pen$a(bt, X1, pmax(pi_t * (1 - pi_t), 1e-8), pi_t)
  mu <- -drop(U %*% solve(M, a_t))                 # Remark 2: mu = U M^{-1} K beta_tilde  ==  -U M^{-1} a(beta_tilde) for the ridge sign convention
  v <- r - mu
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0); Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  Cmat <- matrix(0, G, n); for (g in seq_len(G)) Cmat[g, idx[[g]]] <- 1
  Astar <- Cmat / sqrt(Vg) - U %*% Fi %*% t(X1)
  w_t <- pmax(pi_t * (1 - pi_t), 1e-8)
  list(S_dec = sum(v^2), S_edge = if (is.null(Z)) NA else cf_qf(v, Z), r = r, mu = mu, v = v, Z = Z, beta = beta, beta_tilde = bt,
       Om_MLE = diag(G) - U %*% Fi %*% t(U), R2 = Astar %*% (w_t * t(Astar)), groups = groups, norm_bhat = sqrt(sum(beta[-1]^2)), norm_btilde = sqrt(sum(bt[-1]^2)))
}

# =====================================================================
# Package item D -- ATTENUATION-ADAPTED EDGE DEGREE (Series C Prop 3 as a design rule)
# A degree-k departure survives grouping on the fitted index with non-centrality rho^{2k}.
# rho is observable (Bellec): rho_hat = a_hat / sqrt(a_hat^2 + sigma_hat^2).  Keep degree k while
# rho_hat^{2k} >= tau (default 0.2), at least degree 1.  E8 (2026-09-07): picks the most powerful
# basis in every cell (EDGE-3 at fixed p; EDGE-2 at kappa >= 0.10) and is conservative, not liberal.
# =====================================================================
cf_adaptive_degree <- function(X, y, pc, tau = 0.2, kmax = 3) {
  b <- cf_bellec_adjustments(X, y, pc); rho <- sqrt(b$a2 / (b$a2 + b$s2))
  k <- max(1L, sum(rho^(2 * seq_len(kmax)) >= tau)); list(degree = k, rho_hat = rho)
}
cf_edge_basis <- function(pc, degree) {
  pbar <- vapply(split(seq_along(pc$groups), pc$groups), function(I) mean(pc$pi[I]), 0.0)
  tryCatch(as.matrix(stats::poly(pbar, degree)), error = function(e) NULL)
}
