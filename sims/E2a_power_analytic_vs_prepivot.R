# =====================================================================
# E2a -- Is "more power" already free at fixed p?
#
# Design A (n = 500, p = 5, lambda = 100), the Series C paper's power grid:
# null, omitted quadratic gamma*(x1^2 - 1), omitted interaction gamma*x1*x2,
# gamma in {0.25, 0.5, 0.75, 1}.  On the SAME simulated dataset, three tests:
#   analytic  : closedform.gof, first-order reference, exact weighted chi-square (deflate = 1)
#   prepivot  : shrink.gof (CRAN ebrahim.gof 2.6.0), B = 399, the Series C procedure
#   classical : MLE fit, uncorrected statistic, exact weighted chi-square on Omega_MLE(MLE)
# Paired, so the difference analytic - prepivot has a proper standard error.
#
# Usage: Rscript E2a_power_analytic_vs_prepivot.R [R] [B_boot] [workers]
# =====================================================================
suppressPackageStartupMessages({library(parallel); library(ebrahim.gof)})
args <- commandArgs(trailingOnly = TRUE)
R       <- if (length(args) >= 1) as.integer(args[1]) else 1000
B_BOOT  <- if (length(args) >= 2) as.integer(args[2]) else 399
WORKERS <- if (length(args) >= 3) as.integer(args[3]) else 20
SEED <- 20260906
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))

n <- 500; p <- 5; lambda <- 100; G <- 10
beta0 <- c(0.5, -0.4, 0.3, -0.3, 0.2)
CELLS <- rbind(data.frame(dep = "null", gamma = 0),
               data.frame(dep = "quadratic", gamma = c(0.25, 0.5, 0.75, 1)),
               data.frame(dep = "interaction", gamma = c(0.25, 0.5, 0.75, 1)))

one_rep <- function(r, dep, gamma) {
  set.seed(SEED + r)                                  # same X, y stream within a cell
  X <- matrix(rnorm(n * p), n, p)
  eta <- drop(X %*% beta0) + switch(dep, null = 0, quadratic = gamma * (X[, 1]^2 - 1), interaction = gamma * X[, 1] * X[, 2])
  y <- rbinom(n, 1, cf_expit(eta))
  # analytic (first-order exact reference)
  a <- closedform.gof(X, y, lambda, G = G)
  # prepivoted (the Series C procedure)
  s <- shrink.gof(X, y, lambda = lambda, G = G, basis = c("edge", "decile"), B = B_BOOT, seed = SEED + 10^6 + r)
  # classical: MLE fit, uncorrected statistic, exact reference on Omega_MLE at the MLE
  pm <- cf_pieces(X, y, lambda = 0, G = G)
  data.frame(rep = r, dep = dep, gamma = gamma,
             an_hl = a$SC.HL$p.value, an_edge = a$SC.EDGE$p.value,
             pp_hl = s$SC.HL$p.value, pp_edge = s$SC.EDGE$p.value,
             cl_hl = cf_pvalue(pm$S_dec_unc, pm$Om_MLE), cl_edge = cf_pvalue(pm$S_edge_unc, pm$Om_MLE, pm$Z))
}

t0 <- Sys.time()
cl <- makeCluster(WORKERS); clusterExport(cl, ls())
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(ebrahim.gof)})))
res <- do.call(rbind, lapply(seq_len(nrow(CELLS)), function(k)
  do.call(rbind, parLapply(cl, seq_len(R), one_rep, dep = CELLS$dep[k], gamma = CELLS$gamma[k]))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E2a_power_perrep.csv"), row.names = FALSE)

rej <- function(p) mean(p < 0.05)
agg <- aggregate(cbind(an_hl, an_edge, pp_hl, pp_edge, cl_hl, cl_edge) ~ dep + gamma, data = res, FUN = rej)
# paired differences with SE (over the 0/1 rejection indicators)
pd <- aggregate(cbind(d_hl = (an_hl < 0.05) - (pp_hl < 0.05), d_edge = (an_edge < 0.05) - (pp_edge < 0.05)) ~ dep + gamma, data = res,
                FUN = function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v))))
agg$diff_hl <- pd$d_hl[, "mean"]; agg$se_hl <- pd$d_hl[, "se"]; agg$diff_edge <- pd$d_edge[, "mean"]; agg$se_edge <- pd$d_edge[, "se"]
agg$R <- R; agg$B_boot <- B_BOOT
write.csv(agg, file.path(OUT, "E2a_power_summary.csv"), row.names = FALSE)

options(width = 200)
cat(sprintf("\n=== E2a | design A | R = %d | prepivot B = %d | %.1f min ===\n", R, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
show <- agg; for (v in c("an_hl", "an_edge", "pp_hl", "pp_edge", "cl_hl", "cl_edge")) show[[v]] <- sprintf("%.3f", show[[v]])
show$diff_hl <- sprintf("%+.3f (%.3f)", agg$diff_hl, agg$se_hl); show$diff_edge <- sprintf("%+.3f (%.3f)", agg$diff_edge, agg$se_edge)
print(show[, c("dep", "gamma", "an_hl", "pp_hl", "cl_hl", "diff_hl", "an_edge", "pp_edge", "cl_edge", "diff_edge")], row.names = FALSE)
cat("\nan = analytic first-order exact reference; pp = prepivoted (Series C); cl = classical test on the MLE fit.",
    "\ndiff = analytic - prepivot, paired on identical datasets, (SE).  gamma = 0 rows are size.\n")
