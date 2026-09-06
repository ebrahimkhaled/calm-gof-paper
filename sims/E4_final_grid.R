# E4 -- the definitive size + power grid of the closed-form test vs prepivoting.
# Cells: A (n500,p5,lam100), B kappa=.10 (p40), B kappa=.25 (p100), B kappa=.40 (p160), all lam=416 for B.
# Departures (E1 dial REL): null; index Hermite cubic REL in {.4,.7}; off-index x1*x2 REL in {.4,.7}.
# Tests on the SAME datasets: closed form with reference REF and grouping GRP (args), R2-fitted as the simple variant,
# and prepivoting (CRAN shrink.gof, B bootstrap refits).  Usage: Rscript E4_final_grid.R [R] [B] [workers] [REF] [GRP]
suppressPackageStartupMessages({library(parallel); library(ebrahim.gof)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; B_BOOT <- if (length(args) >= 2) as.integer(args[2]) else 399
W <- if (length(args) >= 3) as.integer(args[3]) else 12; REF <- if (length(args) >= 4) args[4] else "R5"; GRP <- if (length(args) >= 5) args[5] else "loo"
SEED <- 20260908
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416), mk("B", 400, 160, 416))
DEPS <- rbind(data.frame(dep = "null", rel = 0), data.frame(dep = "index", rel = c(.4, .7)), data.frame(dep = "offindex", rel = c(.4, .7)))
one <- function(r, d, dep, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, null = 0, index = d$sd_eta * rel * (Zs^3 - 3 * Zs) / sqrt(6), offindex = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta))
  cf <- closedform.gof(X, y, d$lambda, reference = REF, grouping = GRP)
  c2 <- closedform.gof(X, y, d$lambda, reference = "R2", grouping = "fitted")
  pp <- shrink.gof(X, y, lambda = d$lambda, G = 10, basis = c("edge", "decile"), B = B_BOOT, seed = SEED + 1e6 + r)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, cf_hl = cf$SC.HL$p.value, cf_edge = cf$SC.EDGE$p.value,
             r2_hl = c2$SC.HL$p.value, r2_edge = c2$SC.EDGE$p.value, pp_hl = pp$SC.HL$p.value, pp_edge = pp$SC.EDGE$p.value)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(ebrahim.gof)})))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) do.call(rbind, parLapply(cl, seq_len(R), one, d = d, dep = DEPS$dep[k], rel = DEPS$rel[k]))))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, sprintf("E4_final_%s_%s_perrep.csv", REF, GRP)), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(cf_hl, cf_edge, r2_hl, r2_edge, pp_hl, pp_edge) ~ cell + dep + rel, data = res, FUN = rej)
pd <- aggregate(cbind(d_hl = (cf_hl < .05) - (pp_hl < .05), d_edge = (cf_edge < .05) - (pp_edge < .05)) ~ cell + dep + rel, data = res, FUN = function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
agg$diff_hl <- sprintf("%+.3f (%.3f)", pd$d_hl[, "m"], pd$d_hl[, "se"]); agg$diff_edge <- sprintf("%+.3f (%.3f)", pd$d_edge[, "m"], pd$d_edge[, "se"])
agg$R <- R; agg$B <- B_BOOT; agg$REF <- REF; agg$GRP <- GRP
write.csv(agg, file.path(OUT, sprintf("E4_final_%s_%s_summary.csv", REF, GRP)), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E4 | closed form = %s + %s grouping | R = %d | B = %d | %.1f min ===\n", REF, GRP, R, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("cf_hl", "cf_edge", "r2_hl", "r2_edge", "pp_hl", "pp_edge")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[order(sh$cell, sh$dep, sh$rel), c("cell", "dep", "rel", "cf_hl", "r2_hl", "pp_hl", "diff_hl", "cf_edge", "r2_edge", "pp_edge", "diff_edge")], row.names = FALSE)
cat("\ncf = closed form (REF, GRP); r2 = closed form R2 with fitted grouping; pp = prepivoted. diff = cf - pp, paired. rel = 0 rows are size.\n")
