source("R/closedform_engine.R"); suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
y <- as.numeric(GlaucomaM$Class == "glaucoma"); X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); X1 <- cbind(1, X); D <- diag(c(0, rep(1, ncol(X)))); lam <- 295.3
e1 <- drop(X1 %*% cf_ridge(X1, y, lam, D))
for (m in c(2, 5, 20)) { em <- drop(X1 %*% cf_ridge(X1, y, m * lam, D)); o <- closedform.gof(X, y, lam, grouping = "heavy", heavy_mult = m, inflate = "none", basis = c("decile", "edge"))
  cat(sprintf("glaucoma  mult %2d: Spearman(eta_lam, eta_mlam) = %.3f | same decile: %.0f%% | SC.HL p = %.3f  SC.EDGE p = %.3f\n", m, cor(e1, em, method = "spearman"), 100 * mean(ceiling(rank(e1) / 19.6) == ceiling(rank(em) / 19.6)), o$SC.HL$p.value, o$SC.EDGE$p.value)) }
# the same diagnostic on a design-B draw for contrast
set.seed(5); p <- 100; n <- 400; Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)); Xb <- matrix(rnorm(n * p), n) %*% Sc; yb <- rbinom(n, 1, plogis(Xb %*% b0)); X1b <- cbind(1, Xb); Db <- diag(c(0, rep(1, p)))
eb1 <- drop(X1b %*% cf_ridge(X1b, yb, 416, Db)); eb5 <- drop(X1b %*% cf_ridge(X1b, yb, 5 * 416, Db))
cat(sprintf("design B draw: Spearman(eta_416, eta_2080) = %.3f | same decile: %.0f%%\n", cor(eb1, eb5, method = "spearman"), 100 * mean(ceiling(rank(eb1) / 40) == ceiling(rank(eb5) / 40))))
