# E4b2 -- the FINAL CALM design (degree rule rho2) on E4b's datasets (SEED + r identical): null + Hermite CUBIC rel .4/.7 + x1x2 rel .4/.7.
# E4b's cf_ad column was computed under the retired rho-rule (degree 1 in many cells); this recomputes CALM so the boundary table
# (departures that do not survive attenuation) shows the shipped design.  Prepivot columns are merged from E4b_final_perrep.csv.
# Usage: Rscript E4b2_calm_final_on_E4b.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 6
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
  y <- rbinom(d$n, 1, cf_expit(eta)); cf <- calm.gof(X, y, d$lambda)
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, deg = cf$SC.EDGE.adaptive$degree, rho_hat = cf$SC.EDGE.adaptive$rho_hat,
             calm_hl = cf$SC.HL$p.value, calm_edge = cf$SC.EDGE.adaptive$p.value, calm_e3 = cf$SC.EDGE$p.value)
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d, DEPS$dep[k], DEPS$rel[k]), error = function(e) NULL)))))))
stopCluster(cl); rownames(res) <- NULL
pp <- read.csv(file.path(OUT, "E4b_final_perrep.csv"))[, c("rep", "cell", "dep", "rel", "pp_hl", "pp_edge")]
m <- merge(res, pp, by = c("rep", "cell", "dep", "rel")); write.csv(m, file.path(OUT, "E4b2_final_perrep.csv"), row.names = FALSE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
agg <- aggregate(cbind(calm_hl, calm_edge, calm_e3, pp_hl, pp_edge) ~ cell + dep + rel, m, rej)
pd <- aggregate(cbind(d_hl = (calm_hl < .05) - (pp_hl < .05), d_ed = (calm_edge < .05) - (pp_edge < .05)) ~ cell + dep + rel, m, function(v) c(m = mean(v), se = sd(v) / sqrt(length(v))))
agg$diff_hl <- sprintf("%+.3f (%.3f)", pd$d_hl[, "m"], pd$d_hl[, "se"]); agg$diff_ed <- sprintf("%+.3f (%.3f)", pd$d_ed[, "m"], pd$d_ed[, "se"])
agg$deg <- aggregate(deg ~ cell + dep + rel, m, function(v) as.integer(names(which.max(table(v)))))$deg; agg$R <- R
write.csv(agg, file.path(OUT, "E4b2_final_summary.csv"), row.names = FALSE)
options(width = 220); cat(sprintf("\n=== E4b2 | FINAL CALM on E4b's data (cubic / small x1x2) vs prepivot | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in c("calm_hl", "calm_edge", "calm_e3", "pp_hl", "pp_edge")) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[order(sh$cell, sh$dep, sh$rel), c("cell", "dep", "rel", "deg", "calm_hl", "pp_hl", "diff_hl", "calm_edge", "calm_e3", "pp_edge", "diff_ed")], row.names = FALSE)
