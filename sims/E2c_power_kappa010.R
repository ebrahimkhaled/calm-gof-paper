# =====================================================================
# E2c -- power at kappa = 0.10, where R2 is already VALID (E3: defl 0.99, size .054/.048).
# Design B recipe with p = 40, n = 400, lambda = 416.  Departures matched by the E1 dial
# REL = sd(departure)/sd(eta): index = Hermite cubic (monotone, same ordering),
# off-index = gamma * x1 * x2.  Tests, paired on identical datasets:
#   R2  closed form, variance at pi(beta_tilde)   |   prepivot  shrink.gof (CRAN), B bootstrap refits
# Usage: Rscript E2c_power_kappa010.R [R] [B_boot] [workers]
# =====================================================================
suppressPackageStartupMessages({library(parallel); library(ebrahim.gof)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500
B_BOOT <- if (length(args) >= 2) as.integer(args[2]) else 399
WORKERS <- if (length(args) >= 3) as.integer(args[3]) else 10
SEED <- 20260907
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))

n <- 400; p <- 40; lambda <- 416; G <- 10
Sc <- chol(0.7^abs(outer(1:p, 1:p, "-")))
b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2))
sd_eta <- sqrt(drop(t(b0) %*% (0.7^abs(outer(1:p, 1:p, "-"))) %*% b0))
sd_x1x2 <- sqrt(1 + 0.7^2)
CELLS <- rbind(data.frame(dep = "null", rel = 0),
               data.frame(dep = "index", rel = c(0.4, 0.7)),
               data.frame(dep = "offindex", rel = c(0.4, 0.7)))

one_rep <- function(r, dep, rel) {
  set.seed(SEED + r)
  X <- matrix(rnorm(n * p), n, p) %*% Sc
  eta <- drop(X %*% b0)
  Zs <- eta / sd_eta
  eta <- eta + switch(dep, null = 0,
                      index = sd_eta * rel * (Zs^3 - 3 * Zs) / sqrt(6),
                      offindex = (rel * sd_eta / sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(n, 1, cf_expit(eta))
  pc <- cf_pieces(X, y, lambda, G = G); rf <- cf_references(X, y, lambda, pc)
  s <- shrink.gof(X, y, lambda = lambda, G = G, basis = c("edge", "decile"), B = B_BOOT, seed = SEED + 1e6 + r)
  data.frame(rep = r, dep = dep, rel = rel,
             r2_hl = cf_pvalue(pc$S_dec, rf$R2), r2_edge = cf_pvalue(pc$S_edge, rf$R2, pc$Z),
             r0_hl = cf_pvalue(pc$S_dec, pc$Om_MLE), r0_edge = cf_pvalue(pc$S_edge, pc$Om_MLE, pc$Z),
             pp_hl = s$SC.HL$p.value, pp_edge = s$SC.EDGE$p.value)
}
t0 <- Sys.time()
cl <- makeCluster(WORKERS); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(ebrahim.gof)})))
res <- do.call(rbind, lapply(seq_len(nrow(CELLS)), function(k) do.call(rbind, parLapply(cl, seq_len(R), one_rep, dep = CELLS$dep[k], rel = CELLS$rel[k]))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E2c_power_kappa010_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < 0.05)
agg <- aggregate(cbind(r0_hl, r0_edge, r2_hl, r2_edge, pp_hl, pp_edge) ~ dep + rel, data = res, FUN = rej)
pd <- aggregate(cbind(d_hl = (r2_hl < .05) - (pp_hl < .05), d_edge = (r2_edge < .05) - (pp_edge < .05)) ~ dep + rel, data = res,
                FUN = function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
agg$diff_hl <- sprintf("%+.3f (%.3f)", pd$d_hl[, "m"], pd$d_hl[, "se"]); agg$diff_edge <- sprintf("%+.3f (%.3f)", pd$d_edge[, "m"], pd$d_edge[, "se"])
agg$R <- R; agg$B <- B_BOOT
write.csv(agg, file.path(OUT, "E2c_power_kappa010_summary.csv"), row.names = FALSE)
options(width = 200)
cat(sprintf("\n=== E2c | kappa = 0.10 (n=400, p=40), lambda = 416 | R = %d | B = %d | %.1f min ===\n", R, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("r0_hl", "r0_edge", "r2_hl", "r2_edge", "pp_hl", "pp_edge")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[, c("dep", "rel", "r0_hl", "r2_hl", "pp_hl", "diff_hl", "r0_edge", "r2_edge", "pp_edge", "diff_edge")], row.names = FALSE)
cat("\nr0 = first-order reference; r2 = variance at pi(beta_tilde) (the closed form); pp = prepivoted. diff = r2 - pp, paired. rel = 0 rows are size.\n")
