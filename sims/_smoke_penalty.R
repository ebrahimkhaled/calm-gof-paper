source("R/closedform_engine.R"); set.seed(11)
n <- 500; p <- 5; b0 <- c(.5, -.4, .3, -.3, .2); X <- matrix(rnorm(n * p), n); y <- rbinom(n, 1, plogis(X %*% b0))
# 1) ridge via the general path must equal the dedicated engine
pcR <- cf_pieces(X, y, 100); pcG <- cf_pieces_penalty(X, y, cf_penalty("ridge", lambda = 100))
cat(sprintf("ridge  : S_dec %.6f vs %.6f | S_edge %.6f vs %.6f | max|beta diff| %.2e | max|R2 diff| %.2e\n", pcR$S_dec, pcG$S_dec, pcR$S_edge, pcG$S_edge, max(abs(pcR$beta - pcG$beta)), max(abs(cf_references(X, y, 100, pcR)$R2 - pcG$R2))))
# 2) Firth: fit should match logistf if available; else report the fit and statistics
pf <- cf_pieces_penalty(X, y, cf_penalty("firth"))
cat(sprintf("firth  : ||beta_hat|| %.3f ||beta_tilde|| %.3f (||beta0|| %.3f) | S_dec %.3f (tr Om_MLE %.2f) | p_R0 %.3f p_R2 %.3f | EDGE p_R2 %.3f\n", pf$norm_bhat, pf$norm_btilde, sqrt(sum(b0^2)), pf$S_dec, sum(diag(pf$Om_MLE)), cf_pvalue(pf$S_dec, pf$Om_MLE), cf_pvalue(pf$S_dec, pf$R2), cf_pvalue(pf$S_edge, pf$R2, pf$Z)))
if (requireNamespace("logistf", quietly = TRUE)) { lf <- logistf::logistf(y ~ X, pl = FALSE); cat(sprintf("firth  : max|beta - logistf| = %.2e\n", max(abs(pf$beta - coef(lf))))) } else cat("firth  : logistf not installed (skipping external check)\n")
# 3) generalized ridge
pg <- cf_pieces_penalty(X, y, cf_penalty("gridge", lambda_vec = c(50, 100, 150, 200, 250)))
cat(sprintf("gridge : S_dec %.3f | p_R2 %.3f\n", pg$S_dec, cf_pvalue(pg$S_dec, pg$R2)))
