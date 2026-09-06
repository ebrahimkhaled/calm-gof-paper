# E6 -- Package item G: size of the closed form (R5, fitted grouping) under NON-GAUSSIAN designs.
# Predictor laws: gauss (reference), t3 (heavy tails, scaled to unit variance), binary (Bernoulli(0.3), standardized),
# each with the design's correlation structure imposed through the Cholesky factor (Gaussian copula-free: X = Z_raw %*% Sc).
# Cells: A (n500,p5,lam100), B kappa .10 (p40), B kappa .25 (p100).  Also reports R2 and prepivot-free R0 for context.
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 400; W <- if (length(args) > 1) as.integer(args[2]) else 4
SEED <- 20260909
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
design <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(design("A", 500, 5, 100), design("B", 400, 40, 416), design("B", 400, 100, 416))
LAWS <- c("gauss", "t3", "binary")
draw <- function(n, p, law) switch(law, gauss = matrix(rnorm(n * p), n), t3 = matrix(rt(n * p, 3) / sqrt(3), n), binary = matrix((rbinom(n * p, 1, .3) - .3) / sqrt(.21), n))
one <- function(r, d, law) {
  set.seed(SEED + r); X <- draw(d$n, d$p, law) %*% d$Sc; y <- rbinom(d$n, 1, plogis(drop(X %*% d$b0)))
  pc <- cf_pieces(X, y, d$lambda); r5 <- cf_reference_R5(X, y, d$lambda, pc); rf <- cf_references(X, y, d$lambda, pc)
  data.frame(rep = r, cell = d$name, law = law, S_dec = pc$S_dec, S_edge = pc$S_edge, df5 = cf_implied_df(r5$R5),
             p0_dec = cf_pvalue(pc$S_dec, pc$Om_MLE), p0_edge = cf_pvalue(pc$S_edge, pc$Om_MLE, pc$Z),
             p2_dec = cf_pvalue(pc$S_dec, rf$R2), p2_edge = cf_pvalue(pc$S_edge, rf$R2, pc$Z),
             p5_dec = cf_pvalue(pc$S_dec, r5$R5), p5_edge = cf_pvalue(pc$S_edge, r5$R5, pc$Z), est_sd = r5$est_signal_sd)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(LAWS, function(l) do.call(rbind, parLapply(cl, seq_len(R), one, d = d, law = l)))))); stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E6_nongaussian_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(p0_dec, p0_edge, p2_dec, p2_edge, p5_dec, p5_edge, est_sd) ~ cell + law, data = res, FUN = function(v) if (all(v <= 1 & v >= 0, na.rm = TRUE) && max(v, na.rm = TRUE) <= 1) rej(v) else mean(v))
agg$est_sd <- aggregate(est_sd ~ cell + law, data = res, FUN = mean)$est_sd; agg$defl5 <- aggregate(S_dec ~ cell + law, data = res, FUN = mean)$S_dec / aggregate(df5 ~ cell + law, data = res, FUN = mean)$df5
agg$R <- R; write.csv(agg, file.path(OUT, "E6_nongaussian_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E6 | non-Gaussian designs | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("p0_dec", "p0_edge", "p2_dec", "p2_edge", "p5_dec", "p5_edge", "defl5")) sh[[v]] <- sprintf("%.3f", sh[[v]]); sh$est_sd <- sprintf("%.2f", sh$est_sd)
print(sh[order(sh$cell, sh$law), c("cell", "law", "p0_dec", "p2_dec", "p5_dec", "p0_edge", "p2_edge", "p5_edge", "defl5", "est_sd")], row.names = FALSE)
cat("\np0/p2/p5 = size at 5% under references R0 / R2 / R5 (fitted grouping).  defl5 = mean S / implied df under R5.  true sd(eta): A 0.79, B 1.48-1.50.\n")
