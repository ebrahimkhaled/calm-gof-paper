# E9b -- EDGE with its OWN kappa law, and the adaptive basis.  R5 reference, fitted grouping, 8 cells, H0.
#   edge3        : EDGE-3 on R5, uninflated                 (E9 baseline)
#   edge3_k      : EDGE-3 on (1 + 1.15 kappa) R5            (EDGE-specific law, sims/_delta_law_edge.R)
#   edge2        : EDGE-2 on R5, uninflated                 (the adaptive choice at kappa >= .10)
#   edge2_k      : EDGE-2 on (1 + 0.11 kappa) R5
#   adaptive     : degree by rho_hat rule (keep k while rho_hat^{2k} >= .2), uninflated
#   adaptive_k   : same, inflated by the law of its chosen degree (deg 3: 1.15k; deg 2: 0.11k; deg 1: 0.11k)
#   hl_k         : decile on (1 + 0.24 kappa) R5            (E9's valid decile test, for the record)
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args)) as.integer(args[1]) else 500; W <- if (length(args) > 1) as.integer(args[2]) else 6
SEED <- 20260905
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
design <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p) } else { Sc <- chol(0.7^abs(outer(1:p, 1:p, "-"))); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_n%d_p%d_lam%g", tag, n, p, lambda), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc) }
CELLS <- list(design("A", 500, 5, 100), design("A", 500, 5, 137), design("A", 500, 5, 200), design("B", 400, 100, 50), design("B", 400, 100, 416), design("B", 400, 100, 1000), design("B", 400, 40, 416), design("B", 400, 160, 416))
GLAW <- c(0.11, 0.11, 1.15)   # inflation coefficient on kappa by EDGE degree 1, 2, 3
one <- function(r, d) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; y <- rbinom(d$n, 1, plogis(drop(X %*% d$b0)))
  pc <- cf_pieces(X, y, d$lambda); Sg <- cf_reference_R5(X, y, d$lambda, pc)$R5; kappa <- d$p / d$n
  ad <- cf_adaptive_degree(X, y, pc); Z <- lapply(1:3, function(k) cf_edge_basis(pc, k))
  S <- vapply(1:3, function(k) cf_qf(pc$v, Z[[k]]), 0.0)
  data.frame(rep = r, cell = d$name, kappa = kappa, rho = ad$rho_hat, deg = ad$degree,
    edge3 = cf_pvalue(S[3], Sg, Z[[3]]), edge3_k = cf_pvalue(S[3], (1 + 1.15 * kappa) * Sg, Z[[3]]),
    edge2 = cf_pvalue(S[2], Sg, Z[[2]]), edge2_k = cf_pvalue(S[2], (1 + 0.11 * kappa) * Sg, Z[[2]]),
    adaptive = cf_pvalue(S[ad$degree], Sg, Z[[ad$degree]]), adaptive_k = cf_pvalue(S[ad$degree], (1 + GLAW[ad$degree] * kappa) * Sg, Z[[ad$degree]]),
    hl_k = cf_pvalue(pc$S_dec, (1 + 0.24 * kappa) * Sg))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, parLapply(cl, seq_len(R), one, d = d)))); stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E9b_edge_inflation_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(edge3, edge3_k, edge2, edge2_k, adaptive, adaptive_k, hl_k) ~ cell, data = res, FUN = rej)
agg$deg_mode <- aggregate(deg ~ cell, res, function(v) as.integer(names(which.max(table(v)))))$deg; agg$rho <- aggregate(rho ~ cell, res, median)$rho; agg$R <- R
write.csv(agg, file.path(OUT, "E9b_edge_inflation_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E9b | EDGE-specific inflation + adaptive basis | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("edge3", "edge3_k", "edge2", "edge2_k", "adaptive", "adaptive_k", "hl_k")) sh[[v]] <- sprintf("%.3f", sh[[v]]); sh$rho <- sprintf("%.2f", sh$rho)
print(sh[, c("cell", "rho", "deg_mode", "edge3", "edge3_k", "edge2", "edge2_k", "adaptive", "adaptive_k", "hl_k")], row.names = FALSE)
cat("\nAll on the R5 reference, fitted grouping.  _k = inflated by the basis's own kappa law.  Target .050 (MC SE .010); conservative acceptable.\n")
