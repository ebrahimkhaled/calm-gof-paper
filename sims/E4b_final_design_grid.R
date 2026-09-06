# E4b -- DEFINITIVE size + power grid of the FINAL closed-form design vs prepivoting.
# Final design = closedform.gof(X, y, lambda) defaults: reference R5, fitted grouping, basis-specific kappa-inflation,
# bases SC.HL (decile) and SC.EDGE.adaptive (degree by rho_hat).  Also reports SC.EDGE (EDGE-3 inflated) for the record.
# Cells: A (n500,p5,lam100), B kappa .10 (p40), .25 (p100), .40 (p160), lam=416.  Departures (E1 dial REL): null;
# index Hermite cubic REL {.4,.7}; off-index x1*x2 REL {.4,.7}.  Paired with prepivoting (CRAN shrink.gof, B refits).
# Usage: Rscript E4b_final_design_grid.R [R] [B] [workers]
suppressPackageStartupMessages({library(parallel); library(ebrahim.gof)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; B_BOOT <- if (length(args) >= 2) as.integer(args[2]) else 399; W <- if (length(args) >= 3) as.integer(args[3]) else 12
SEED <- 20260908                       # SAME stream as E4 -> identical datasets
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
  cf <- closedform.gof(X, y, d$lambda)                      # the final design, all defaults
  pp <- shrink.gof(X, y, lambda = d$lambda, G = 10, basis = c("edge", "decile"), B = B_BOOT, seed = SEED + 1e6 + r)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, deg = cf$SC.EDGE.adaptive$degree,
             cf_hl = cf$SC.HL$p.value, cf_ad = cf$SC.EDGE.adaptive$p.value, cf_e3 = cf$SC.EDGE$p.value,
             pp_hl = pp$SC.HL$p.value, pp_edge = pp$SC.EDGE$p.value)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(ebrahim.gof)})))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) do.call(rbind, parLapply(cl, seq_len(R), one, d = d, dep = DEPS$dep[k], rel = DEPS$rel[k]))))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E4b_final_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(cf_hl, cf_ad, cf_e3, pp_hl, pp_edge) ~ cell + dep + rel, data = res, FUN = rej)
pd <- aggregate(cbind(d_hl = (cf_hl < .05) - (pp_hl < .05), d_ad = (cf_ad < .05) - (pp_edge < .05)) ~ cell + dep + rel, data = res, FUN = function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
agg$diff_hl <- sprintf("%+.3f (%.3f)", pd$d_hl[, "m"], pd$d_hl[, "se"]); agg$diff_ad <- sprintf("%+.3f (%.3f)", pd$d_ad[, "m"], pd$d_ad[, "se"])
agg$deg <- aggregate(deg ~ cell + dep + rel, data = res, FUN = function(v) as.integer(names(which.max(table(v)))))$deg
agg$R <- R; agg$B <- B_BOOT; write.csv(agg, file.path(OUT, "E4b_final_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E4b | FINAL closed form vs prepivoting | R = %d | B = %d | %.1f min ===\n", R, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("cf_hl", "cf_ad", "cf_e3", "pp_hl", "pp_edge")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[order(sh$cell, sh$dep, sh$rel), c("cell", "dep", "rel", "deg", "cf_hl", "pp_hl", "diff_hl", "cf_ad", "cf_e3", "pp_edge", "diff_ad")], row.names = FALSE)
cat("\ncf_hl = SC.HL (kappa-inflated R5); cf_ad = SC.EDGE.adaptive; cf_e3 = SC.EDGE-3 inflated; pp = prepivoted.  diff = closed form - prepivot, paired (SE).  rel = 0 rows are size.\n")
