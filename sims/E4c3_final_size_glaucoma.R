# E4c3 -- SIZE of the REDESIGNED CALM (degree_rule = "kappa": EDGE-3 at fixed p, EDGE-2 in the proportional regime; HL kappa-inflated R5)
# in the eight Table-1 cells, plus the glaucoma p-values and timing under the new default.  No bootstrap.  R = 500 (Table 1 convention).
# Usage: Rscript E4c3_final_size_glaucoma.R [R] [workers]
suppressPackageStartupMessages({library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260908
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d_lam%d", tag, p, lambda), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(mk("A", 500, 5, 100), mk("A", 500, 5, 137), mk("A", 500, 5, 200), mk("B", 400, 40, 416), mk("B", 400, 100, 50), mk("B", 400, 100, 416), mk("B", 400, 100, 1000), mk("B", 400, 160, 416))
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, cf_expit(drop(X %*% d$b0)))
  cf <- calm.gof(X, y, d$lambda)
  data.frame(rep = r, cell = d$name, hl = cf$SC.HL$p.value, edge_k = cf$SC.EDGE.adaptive$p.value, deg = cf$SC.EDGE.adaptive$degree, rho_hat = cf$SC.EDGE.adaptive$rho_hat, e3 = cf$SC.EDGE$p.value)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one, d = d))))
stopCluster(cl); write.csv(res, file.path(OUT, "E4c3_size_perrep.csv"), row.names = FALSE)
agg <- aggregate(cbind(hl, edge_k, e3) ~ cell, res, function(v) mean(v < .05)); agg$deg <- aggregate(deg ~ cell, res, function(v) v[1])$deg; agg$rho_hat <- round(aggregate(rho_hat ~ cell, res, mean)$rho_hat, 2); agg$R <- R
write.csv(agg, file.path(OUT, "E4c3_size_summary.csv"), row.names = FALSE)
cat(sprintf("\n=== E4c3 | size of the redesigned CALM | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins")))); print(agg, digits = 3, row.names = FALSE)
# ---- glaucoma at the cross-validated penalty (Series C: lambda = 295, G = 10) ----
suppressPackageStartupMessages(library(TH.data)); data("GlaucomaM", package = "TH.data")
Xg <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); yg <- as.numeric(GlaucomaM$Class == "glaucoma")   # standardised, as glmnet and Series C
tt <- system.time(cg <- calm.gof(Xg, yg, 295.3))
cat(sprintf("\nGlaucoma (n = %d, p = %d, kappa = %.2f, lambda = 295): SC.HL p = %.4f | SC.EDGE degree-%d p = %.4f (rho_hat %.2f) | SC.EDGE-3 p = %.4f | %.2f s\n",
            nrow(Xg), ncol(Xg), ncol(Xg) / nrow(Xg), cg$SC.HL$p.value, cg$SC.EDGE.adaptive$degree, cg$SC.EDGE.adaptive$p.value, cg$SC.EDGE.adaptive$rho_hat, cg$SC.EDGE$p.value, tt[["elapsed"]]))
writeLines(sprintf("%s,%g", c("hl", "edge_k", "degree", "rho_hat", "e3", "seconds"), c(cg$SC.HL$p.value, cg$SC.EDGE.adaptive$p.value, cg$SC.EDGE.adaptive$degree, cg$SC.EDGE.adaptive$rho_hat, cg$SC.EDGE$p.value, tt[["elapsed"]])), file.path(OUT, "E4c3_glaucoma.csv"))
