# E7 -- Package item F: the closed form under OTHER smooth penalties (Series C Remark 2, made numerical).
# Penalties: firth (Jeffreys score; the clinical standard for separation), gridge (generalized ridge with
# lambda_j = lambda * (1 + j/p), i.e. unequal weights), and ridge as the control.
# Per cell and penalty, under H0: the UNCORRECTED statistic referred to Omega_MLE (does the penalty distort the
# textbook test?), and the corrected statistic under R0 (first order) and R2 (variance at the debiased generator).
# Designs: A (n500,p5,lam100) and B (n400,p100,lam416).  Firth has no lambda.
# Usage: Rscript E7_penalty_generality.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 400; W <- if (length(args) > 1) as.integer(args[2]) else 10
SEED <- 20260910
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
design <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(design("A", 500, 5, 100), design("B", 400, 100, 416))
PENS <- c("ridge", "gridge", "firth")
one <- function(r, d, pen_name) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, plogis(drop(X %*% d$b0)))
  pen <- switch(pen_name, ridge = cf_penalty("ridge", lambda = d$lambda), gridge = cf_penalty("gridge", lambda_vec = d$lambda * (1 + seq_len(d$p) / d$p)), firth = cf_penalty("firth"))
  pc <- tryCatch(cf_pieces_penalty(X, y, pen), error = function(e) NULL)
  if (is.null(pc)) return(data.frame(rep = r, cell = d$name, pen = pen_name, S_unc = NA, S_dec = NA, S_edge = NA, tr = NA, p_unc = NA, p_R0 = NA, p_R2 = NA, pe_unc = NA, pe_R0 = NA, pe_R2 = NA, nb = NA, nbt = NA))
  data.frame(rep = r, cell = d$name, pen = pen_name, S_unc = sum(pc$r^2), S_dec = pc$S_dec, S_edge = pc$S_edge, tr = sum(diag(pc$Om_MLE)),
             p_unc = cf_pvalue(sum(pc$r^2), pc$Om_MLE), p_R0 = cf_pvalue(pc$S_dec, pc$Om_MLE), p_R2 = cf_pvalue(pc$S_dec, pc$R2),
             pe_unc = cf_pvalue(cf_qf(pc$r, pc$Z), pc$Om_MLE, pc$Z), pe_R0 = cf_pvalue(pc$S_edge, pc$Om_MLE, pc$Z), pe_R2 = cf_pvalue(pc$S_edge, pc$R2, pc$Z),
             nb = pc$norm_bhat, nbt = pc$norm_btilde)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(PENS, function(pn) do.call(rbind, parLapply(cl, seq_len(R), one, d = d, pen_name = pn)))))); stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E7_penalty_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- do.call(rbind, lapply(split(res, list(res$cell, res$pen), drop = TRUE), function(s) data.frame(cell = s$cell[1], pen = s$pen[1], fails = sum(is.na(s$S_dec)),
  size_unc = rej(s$p_unc), size_R0 = rej(s$p_R0), size_R2 = rej(s$p_R2), esize_unc = rej(s$pe_unc), esize_R0 = rej(s$pe_R0), esize_R2 = rej(s$pe_R2),
  defl_R0 = mean(s$S_dec, na.rm = TRUE) / mean(s$tr, na.rm = TRUE), norm_bhat = mean(s$nb, na.rm = TRUE), norm_btilde = mean(s$nbt, na.rm = TRUE))))
agg$R <- R; write.csv(agg, file.path(OUT, "E7_penalty_summary.csv"), row.names = FALSE)
options(width = 200); cat(sprintf("\n=== E7 | penalty generality | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("size_unc", "size_R0", "size_R2", "esize_unc", "esize_R0", "esize_R2", "defl_R0")) sh[[v]] <- sprintf("%.3f", sh[[v]]); for (v in c("norm_bhat", "norm_btilde")) sh[[v]] <- sprintf("%.2f", sh[[v]])
print(sh[order(sh$cell, sh$pen), ], row.names = FALSE)
cat("\nsize_unc = UNCORRECTED statistic on Omega_MLE (the textbook test under this penalty); size_R0 / size_R2 = corrected, first-order / debiased-variance reference; esize_* = EDGE basis. ||beta0||: A 0.79, B 2.31.\n")
