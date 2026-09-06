# E4c -- POWER of the final design vs prepivoting against departures that SURVIVE attenuation.
# E4b used a Hermite CUBIC at rel .4/.7: at kappa the adaptive basis (k = 1-2) is orthogonal to a cubic by construction and
# rho^6 ~ .1 kills the rest, so E4b's kappa rows were a null-power design (size table only).  Here: Hermite QUADRATIC along the
# index at rel {.5, 1, 1.5} (E8 showed EDGE-2 has .5-.7 power against it at kappa .25) and the off-index x1*x2 at rel {1, 1.5}.
# Same cells, same SEED stream as E4b (rep r -> identical X, eta), paired with prepivoting on identical data.
# Usage: Rscript E4c_power_survivors.R [R] [B] [workers]
suppressPackageStartupMessages({library(parallel); library(ebrahim.gof)})
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 400; B_BOOT <- if (length(args) >= 2) as.integer(args[2]) else 299; W <- if (length(args) >= 3) as.integer(args[3]) else 8
SEED <- 20260908
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416), mk("B", 400, 160, 416))
DEPS <- rbind(data.frame(dep = "quadratic", rel = c(.5, 1, 1.5)), data.frame(dep = "offindex", rel = c(1, 1.5)))
one <- function(r, d, dep, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, quadratic = d$sd_eta * rel * (Zs^2 - 1) / sqrt(2), offindex = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta))
  cf <- calm.gof(X, y, d$lambda)
  pp <- shrink.gof(X, y, lambda = d$lambda, G = 10, basis = c("edge", "decile"), B = B_BOOT, seed = SEED + 1e6 + r)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, deg = cf$SC.EDGE.adaptive$degree, rho_hat = cf$SC.EDGE.adaptive$rho_hat,
             cf_hl = cf$SC.HL$p.value, cf_ad = cf$SC.EDGE.adaptive$p.value, cf_e3 = cf$SC.EDGE$p.value,
             pp_hl = pp$SC.HL$p.value, pp_edge = pp$SC.EDGE$p.value)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(CompQuadForm); library(ebrahim.gof)})))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) {
  out <- do.call(rbind, parLapply(cl, seq_len(R), one, d = d, dep = DEPS$dep[k], rel = DEPS$rel[k]))
  write.csv(out, file.path(OUT, sprintf("E4c_perrep_%s_%s_%g.csv", d$name, DEPS$dep[k], DEPS$rel[k])), row.names = FALSE)   # survives a crash
  cat(sprintf("done %s %s rel %.1f  (%.0f min)\n", d$name, DEPS$dep[k], DEPS$rel[k], as.numeric(difftime(Sys.time(), t0, units = "mins")))); out }))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E4c_power_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(cf_hl, cf_ad, cf_e3, pp_hl, pp_edge) ~ cell + dep + rel, data = res, FUN = rej)
pd <- aggregate(cbind(d_hl = (cf_hl < .05) - (pp_hl < .05), d_ad = (cf_ad < .05) - (pp_edge < .05)) ~ cell + dep + rel, data = res, FUN = function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
agg$diff_hl <- sprintf("%+.3f (%.3f)", pd$d_hl[, "m"], pd$d_hl[, "se"]); agg$diff_ad <- sprintf("%+.3f (%.3f)", pd$d_ad[, "m"], pd$d_ad[, "se"])
agg$deg <- aggregate(deg ~ cell + dep + rel, data = res, FUN = function(v) as.integer(names(which.max(table(v)))))$deg
agg$rho_hat <- round(aggregate(rho_hat ~ cell + dep + rel, data = res, FUN = mean)$rho_hat, 2)
agg$R <- R; agg$B <- B_BOOT; write.csv(agg, file.path(OUT, "E4c_power_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E4c | power vs prepivoting, surviving departures | R = %d | B = %d | %.1f min ===\n", R, B_BOOT, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("cf_hl", "cf_ad", "cf_e3", "pp_hl", "pp_edge")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[order(sh$cell, sh$dep, sh$rel), c("cell", "dep", "rel", "deg", "rho_hat", "cf_hl", "pp_hl", "diff_hl", "cf_ad", "cf_e3", "pp_edge", "diff_ad")], row.names = FALSE)
cat("\ncf_hl = SC.HL (CALM); cf_ad = SC.EDGE.adaptive (CALM); cf_e3 = SC.EDGE-3 inflated; pp = prepivoted.  diff = CALM - prepivot, paired (SE).\n")
