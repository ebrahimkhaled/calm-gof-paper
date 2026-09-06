# Glaucoma diagnostics for the referee point "CALM p = .012 < Sigma_d p = .024 despite inflation => Sigma_5 < Sigma_d here: explain".
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "..", "R", "closedform_engine.R"))
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); y <- as.numeric(GlaucomaM$Class == "glaucoma"); lam <- 295.3
pc <- cf_pieces(X, y, lam, G = 10); refs <- cf_references(X, y, lam, pc); r5 <- cf_reference_R5(X, y, lam, pc); b <- cf_bellec_adjustments(X, y, pc)
Sd <- refs$R2; S5 <- r5$R5; S0 <- pc$Om_MLE
X1 <- cbind(1, X); D <- diag(c(0, rep(1, ncol(X)))); bh <- cf_ridge(X1, y, lam, D); eh <- drop(X1 %*% bh)
Fm <- crossprod(X1, (cf_expit(eh) * (1 - cf_expit(eh))) * X1); bt <- bh + solve(Fm, lam * D %*% bh); et <- drop(X1 %*% bt)
cat(sprintf("n=%d p=%d kappa=%.2f lambda=%.1f\n", nrow(X), ncol(X), ncol(X) / nrow(X), lam))
cat(sprintf("sd(eta_hat)=%.3f  sd(eta_tilde)=%.3f  |beta_hat|=%.2f |beta_tilde|=%.2f\n", sd(eh), sd(et), sqrt(sum(bh[-1]^2)), sqrt(sum(bt[-1]^2))))
cat(sprintf("Bellec: a_hat=%.3f sigma_hat=%.3f rho_hat=%.3f\n", sqrt(b$a2), sqrt(b$s2), sqrt(b$a2 / (b$a2 + b$s2))))
str(r5[setdiff(names(r5), "R5")], give.attr = FALSE)
cat(sprintf("tr Omega_MLE=%.3f  tr Sigma_d=%.3f  tr Sigma_5=%.3f  ratio S5/Sd=%.3f\n", sum(diag(S0)), sum(diag(Sd)), sum(diag(S5)), sum(diag(S5)) / sum(diag(Sd))))
w0 <- r5$w0 %||% NULL
wt <- cf_expit(et) * (1 - cf_expit(et)); wh <- cf_expit(eh) * (1 - cf_expit(eh))
cat(sprintf("mean w_hat=%.4f  mean w_tilde=%.4f", mean(wh), mean(wt)))
if (!is.null(r5$w0)) cat(sprintf("  mean w0_hat=%.4f", mean(r5$w0))); cat("\n")
kap <- ncol(X) / nrow(X)
cat(sprintf("p-values SC.HL: Omega_MLE %.4f | Sigma_d %.4f | Sigma_5 %.4f | CALM (1+.24k)Sigma_5 %.4f\n", cf_pvalue(pc$S_dec, S0), cf_pvalue(pc$S_dec, Sd), cf_pvalue(pc$S_dec, S5), cf_pvalue(pc$S_dec, (1 + .24 * kap) * S5)))
cat(sprintf("p-values SC.EDGE-3: Omega_MLE %.4f | Sigma_d %.4f | Sigma_5 %.4f | CALM (1+1.15k)Sigma_5 %.4f\n", cf_pvalue(pc$S_edge, S0, pc$Z), cf_pvalue(pc$S_edge, Sd, pc$Z), cf_pvalue(pc$S_edge, S5, pc$Z), cf_pvalue(pc$S_edge, (1 + 1.15 * kap) * S5, pc$Z)))
Z2 <- cf_edge_basis(pc, 2); cat(sprintf("EDGE-2 under CALM: %.4f ; EDGE-1 (slope direction): %.4f\n", cf_pvalue(cf_qf(pc$v, Z2), (1 + .11 * kap) * S5, Z2), cf_pvalue(cf_qf(pc$v, cf_edge_basis(pc, 1)), (1 + .11 * kap) * S5, cf_edge_basis(pc, 1))))
cat("per-decile corrected residual v:", round(pc$v, 2), "\n")
cat("EDGE-3 components (coef on degree 1,2,3):", round(drop(crossprod(pc$Z, pc$v)), 2), "\n")
