# E4c2 -- CALM design variants on the SAME data as E4b/E4c (SEED + r regenerates X, eta, y exactly), no bootstrap, fast.
# Why: E4c interim showed the adaptive rule picking degree 1 UNDER ALTERNATIVES at kappa (rho_hat collapses when the model
# is misspecified, so the rule discards the quadratic direction that carries the departure) -> adaptive EDGE power .17-.44
# vs prepivot EDGE-3 .46-.88 at kappa .10/.25.  Also HL(CALM) < HL(prepivot) by ~10 pp at kappa .25.
# Variants computed per rep (all p-values):
#   hl_R5_inf   = SC.HL, R5, kappa-inflation (CALM as shipped)      hl_R5     = SC.HL, R5, no inflation
#   hl_R2       = SC.HL, variance at beta-tilde, no inflation       hl_R0     = first order
#   e2_R5_inf   = EDGE-2 fixed, R5, e2-inflation (0.11 kappa)       e2_R5     = EDGE-2, R5, no inflation
#   e2_R2       = EDGE-2, R2                                         e3_R5_inf = EDGE-3, R5, e3-inflation (CALM SC.EDGE)
#   e3_R2       = EDGE-3, R2, no inflation                           ad_R5_inf = adaptive (as shipped), + its degree and rho_hat
#   kd_R5_inf   = "degree from kappa" rule: k = 3 if p/n < .05 else 2, R5, basis-specific inflation
# Cells: A, kappa .10/.25/.40; deps: null (size), quadratic rel .5/1/1.5, offindex rel 1/1.5.  Usage: Rscript E4c2_calm_variants.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 400; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260908
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mk <- function(tag, n, p, lambda) { if (tag == "A") { b0 <- c(.5, -.4, .3, -.3, .2); Sc <- diag(p); Sig <- diag(p) } else { Sig <- 0.7^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p); b0 <- b0 * 2.31 / sqrt(sum(b0^2)) }
  list(name = sprintf("%s_p%d", tag, p), n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, sd_eta = sqrt(drop(t(b0) %*% Sig %*% b0)), sd_x1x2 = if (tag == "A") 1 else sqrt(1 + .49)) }
CELLS <- list(mk("A", 500, 5, 100), mk("B", 400, 40, 416), mk("B", 400, 100, 416), mk("B", 400, 160, 416))
DEPS <- rbind(data.frame(dep = "null", rel = 0), data.frame(dep = "quadratic", rel = c(.5, 1, 1.5)), data.frame(dep = "offindex", rel = c(1, 1.5)))
one <- function(r, d, dep, rel) {
  set.seed(SEED + r); X <- matrix(rnorm(d$n * d$p), d$n) %*% d$Sc; eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + switch(dep, null = 0, quadratic = d$sd_eta * rel * (Zs^2 - 1) / sqrt(2), offindex = (rel * d$sd_eta / d$sd_x1x2) * X[, 1] * X[, 2])
  y <- rbinom(d$n, 1, cf_expit(eta)); n <- d$n; kappa <- d$p / n
  pc <- cf_pieces(X, y, d$lambda, G = 10)
  R5 <- cf_reference_R5(X, y, d$lambda, pc)$R5; R2 <- cf_references(X, y, d$lambda, pc)$R2; R0 <- pc$Om_MLE
  g <- function(deg) kappa * c(dec = 0.24, e1 = 0.11, e2 = 0.11, e3 = 1.15)[[deg]]
  Z2 <- cf_edge_basis(pc, 2); S2 <- cf_qf(pc$v, Z2); Z3 <- pc$Z; S3 <- pc$S_edge
  ad <- cf_adaptive_degree(X, y, pc); Zk <- cf_edge_basis(pc, ad$degree); Sk <- cf_qf(pc$v, Zk); gk <- g(c("e1", "e2", "e3")[ad$degree])
  kd <- if (kappa < .05) 3L else 2L; Zkd <- if (kd == 3L) Z3 else Z2; Skd <- if (kd == 3L) S3 else S2; gkd <- g(c("e1", "e2", "e3")[kd])
  data.frame(rep = r, cell = d$name, dep = dep, rel = rel, deg_ad = ad$degree, rho_hat = ad$rho_hat,
             hl_R5_inf = cf_pvalue(pc$S_dec, (1 + g("dec")) * R5), hl_R5 = cf_pvalue(pc$S_dec, R5), hl_R2 = cf_pvalue(pc$S_dec, R2), hl_R0 = cf_pvalue(pc$S_dec, R0),
             e2_R5_inf = cf_pvalue(S2, (1 + g("e2")) * R5, Z2), e2_R5 = cf_pvalue(S2, R5, Z2), e2_R2 = cf_pvalue(S2, R2, Z2),
             e3_R5_inf = cf_pvalue(S3, (1 + g("e3")) * R5, Z3), e3_R2 = cf_pvalue(S3, R2, Z3),
             ad_R5_inf = cf_pvalue(Sk, (1 + gk) * R5, Zk), kd_R5_inf = cf_pvalue(Skd, (1 + gkd) * R5, Zkd))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
res <- do.call(rbind, lapply(CELLS, function(d) do.call(rbind, lapply(seq_len(nrow(DEPS)), function(k) {
  out <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one(r, d, DEPS$dep[k], DEPS$rel[k]), error = function(e) NULL)))
  cat(sprintf("done %s %s rel %.1f  n=%d  (%.1f min)\n", d$name, DEPS$dep[k], DEPS$rel[k], nrow(out), as.numeric(difftime(Sys.time(), t0, units = "mins")))); out }))))
stopCluster(cl); rownames(res) <- NULL
write.csv(res, file.path(OUT, "E4c2_variants_perrep.csv"), row.names = FALSE)
# merge the prepivot columns from E4b (null) and E4c (departures): identical data by construction
pp <- rbind(read.csv(file.path(OUT, "E4b_final_perrep.csv"))[, c("rep", "cell", "dep", "rel", "pp_hl", "pp_edge")],
            do.call(rbind, lapply(list.files(OUT, "^E4c_perrep_", full.names = TRUE), function(f) read.csv(f)[, c("rep", "cell", "dep", "rel", "pp_hl", "pp_edge")])))
m <- merge(res, pp, by = c("rep", "cell", "dep", "rel"), all.x = TRUE)
rej <- function(v) mean(v < .05, na.rm = TRUE)
vars <- c("hl_R5_inf", "hl_R5", "hl_R2", "hl_R0", "pp_hl", "e2_R5_inf", "e2_R5", "e2_R2", "e3_R5_inf", "e3_R2", "ad_R5_inf", "kd_R5_inf", "pp_edge")
agg <- aggregate(m[, vars], by = m[, c("cell", "dep", "rel")], FUN = rej)
agg$deg_ad <- aggregate(deg_ad ~ cell + dep + rel, m, function(v) as.integer(names(which.max(table(v)))))$deg_ad
agg$rho_hat <- round(aggregate(rho_hat ~ cell + dep + rel, m, mean)$rho_hat, 2); agg$R <- R
write.csv(agg, file.path(OUT, "E4c2_variants_summary.csv"), row.names = FALSE)
options(width = 250); cat(sprintf("\n=== E4c2 | CALM variants vs prepivot on identical data | R = %d | %.1f min ===\n", R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
sh <- agg; for (v in vars) sh[[v]] <- sprintf("%.3f", sh[[v]])
print(sh[order(sh$cell, sh$dep, sh$rel), c("cell", "dep", "rel", "deg_ad", "rho_hat", "hl_R5_inf", "hl_R5", "hl_R2", "pp_hl", "e2_R5_inf", "e2_R5", "e2_R2", "e3_R5_inf", "ad_R5_inf", "kd_R5_inf", "pp_edge")], row.names = FALSE)
cat("\nrel = 0 rows are SIZE.  _inf = kappa-inflated.  kd = degree-from-kappa rule (3 at fixed p, 2 in the proportional regime).\n")
