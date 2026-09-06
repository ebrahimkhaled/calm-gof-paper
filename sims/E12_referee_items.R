# E12 -- the simulation items the editor's report asks for, under the FINAL degree rule "rho2k"
#   (A) size of CALM in the 8 calibration cells + the 17 held-out cells, under the new rule
#   (B) power by basis against cubic and quadratic departures at rel 1 and 1.5 (does rho2k fix the E8 hole at fixed p?)
#   (C) tau sensitivity: tau in {0.05, 0.1, 0.2, 0.3} -> size and chosen degree
#   (D) oracle-Sigma_x check: recompute the Bellec adjustments with the TRUE Sigma_x at kappa .25 and .40
#   (E) scale recovery per cell: estimated vs true sd of the linear predictor
#   (F) uninflated EDGE-2 size (is the EDGE inflation doing anything?)
# Usage: Rscript E12_referee_items.R [R] [workers]
suppressPackageStartupMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args) >= 1) as.integer(args[1]) else 500; W <- if (length(args) >= 2) as.integer(args[2]) else 8
SEED <- 20260910
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))); ROOT <- dirname(HERE); OUT <- file.path(ROOT, "results")
source(file.path(ROOT, "R", "closedform_engine.R"))
mkB <- function(name, n, p, lambda, ar = 0.7, norm = 2.31, law = "gaussian", G = 10) {
  Sig <- if (ar == 0) diag(p) else ar^abs(outer(1:p, 1:p, "-")); Sc <- chol(Sig); b0 <- rep_len(c(.35, -.3, .25, -.2, .15), p)
  b0 <- b0 * norm / sqrt(drop(t(b0) %*% Sig %*% b0))
  list(name = name, n = n, p = p, lambda = lambda, b0 = b0, Sc = Sc, Sig = Sig, law = law, G = G, sd_eta = norm) }
mkA <- function(name, lambda = 100, law = "gaussian") { b <- c(.5, -.4, .3, -.3, .2)
  list(name = name, n = 500, p = 5, lambda = lambda, b0 = b, Sc = diag(5), Sig = diag(5), law = law, G = 10, sd_eta = sqrt(sum(b^2))) }
CAL <- list(mkA("A_lam100", 100), mkA("A_lam137", 137), mkA("A_lam200", 200), mkB("B_k10", 400, 40, 416),
            mkB("B_k25_lam50", 400, 100, 50), mkB("B_k25_lam416", 400, 100, 416), mkB("B_k25_lam1000", 400, 100, 1000), mkB("B_k40", 400, 160, 416))
HELD <- list(mkB("H1_k15", 400, 60, 416), mkB("H2_k33", 400, 132, 416), mkB("H3_k40_lam50", 400, 160, 50), mkB("H4_k40_lam1000", 400, 160, 1000),
             mkB("H5_k25_ar03", 400, 100, 416, ar = 0.3), mkB("H6_k25_weak", 400, 100, 416, norm = 1.5), mkB("H7_k25_strong", 400, 100, 416, norm = 3.0),
             mkB("H8_n800_k25", 800, 200, 832), mkB("H9_k25_G20", 400, 100, 416, G = 20),
             mkA("A_t3", 100, "t3"), mkA("A_binary", 100, "binary"), mkB("k10_t3", 400, 40, 416, law = "t3"), mkB("k10_binary", 400, 40, 416, law = "binary"),
             mkB("k25_t3", 400, 100, 416, law = "t3"), mkB("k25_binary", 400, 100, 416, law = "binary"), mkB("k40_t3", 400, 160, 416, law = "t3"), mkB("k40_binary", 400, 160, 416, law = "binary"))
draw <- function(d, r) { set.seed(SEED + r)
  Z <- switch(d$law, gaussian = matrix(rnorm(d$n * d$p), d$n), t3 = matrix(rt(d$n * d$p, 3) / sqrt(3), d$n), binary = matrix((rbinom(d$n * d$p, 1, .3) - .3) / sqrt(.21), d$n))
  Z %*% d$Sc }
TAUS <- c(0.05, 0.10, 0.20, 0.30)
one_size <- function(r, d) {
  X <- draw(d, r); eta0 <- drop(X %*% d$b0); y <- rbinom(d$n, 1, cf_expit(eta0))
  pc <- cf_pieces(X, y, d$lambda, G = d$G); r5 <- cf_reference_R5(X, y, d$lambda, pc); R5 <- r5$R5
  ad <- cf_adaptive_degree(X, y, pc); kappa <- d$p / d$n
  g <- function(deg) kappa * c(dec = 0.24, e1 = 0.11, e2 = 0.11, e3 = 1.15)[[deg]]
  pk <- function(k) { Zk <- cf_edge_basis(pc, k); cf_pvalue(cf_qf(pc$v, Zk), (1 + g(c("e1","e2","e3")[k])) * R5, Zk) }
  p2 <- pk(2); p3 <- pk(3)
  Z2 <- cf_edge_basis(pc, 2); p2_noinf <- cf_pvalue(cf_qf(pc$v, Z2), R5, Z2)                     # (F)
  taudeg <- sapply(TAUS, function(tt) { k <- max(1L, sum(ad$rho_hat^(2 * 1:3) >= tt)); if (kappa < 0.05 || k >= 3L) 3L else 2L })  # (C) under rho2k
  out <- data.frame(rep = r, cell = d$name, kappa = kappa, rho_hat = ad$rho_hat, deg_final = if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L,
                    hl = cf_pvalue(pc$S_dec, (1 + g("dec")) * R5), e2 = p2, e3 = p3, e2_noinf = p2_noinf,
                    est_sd = r5$est_signal_sd, true_sd = d$sd_eta)                                # (E)
  for (i in seq_along(TAUS)) { out[[paste0("tau", i)]] <- if (taudeg[i] == 3L) p3 else p2; out[[paste0("tdeg", i)]] <- taudeg[i] }
  out$edge <- if (out$deg_final == 3L) p3 else p2
  out
}
one_pow <- function(r, d, deg, rel) {
  X <- draw(d, r); eta <- drop(X %*% d$b0); Zs <- eta / d$sd_eta
  eta <- eta + d$sd_eta * rel * switch(as.character(deg), "2" = (Zs^2 - 1) / sqrt(2), "3" = (Zs^3 - 3 * Zs) / sqrt(6))
  y <- rbinom(d$n, 1, cf_expit(eta))
  pc <- cf_pieces(X, y, d$lambda, G = d$G); R5 <- cf_reference_R5(X, y, d$lambda, pc)$R5
  ad <- cf_adaptive_degree(X, y, pc); kappa <- d$p / d$n
  g <- function(dg) kappa * c(dec = 0.24, e1 = 0.11, e2 = 0.11, e3 = 1.15)[[dg]]
  pk <- function(k) { Zk <- cf_edge_basis(pc, k); cf_pvalue(cf_qf(pc$v, Zk), (1 + g(c("e1","e2","e3")[k])) * R5, Zk) }
  kf <- if (kappa < 0.05 || ad$degree >= 3L) 3L else 2L
  data.frame(rep = r, cell = d$name, dep_deg = deg, rel = rel, rho_hat = ad$rho_hat, deg_final = kf,
             hl = cf_pvalue(pc$S_dec, (1 + g("dec")) * R5), e2 = pk(2), e3 = pk(3), edge = if (kf == 3L) pk(3) else pk(2),
             edge_old = { ko <- max(2L, ad$degree); if (ko == 3L) pk(3) else pk(2) })
}
one_oracle <- function(r, d) {                                                                   # (D)
  X <- draw(d, r); y <- rbinom(d$n, 1, cf_expit(drop(X %*% d$b0)))
  pc <- cf_pieces(X, y, d$lambda, G = d$G)
  b_hat <- cf_bellec_adjustments(X, y, pc)
  b_or <- cf_bellec_adjustments(X, y, pc, Sigma_x = d$Sig)
  data.frame(rep = r, cell = d$name, a_hat = sqrt(b_hat$a2), a_or = sqrt(b_or$a2), s_hat = sqrt(b_hat$s2), s_or = sqrt(b_or$s2),
             rho_hat = sqrt(b_hat$a2 / (b_hat$a2 + b_hat$s2)), rho_or = sqrt(b_or$a2 / (b_or$a2 + b_or$s2)))
}
t0 <- Sys.time(); cl <- makeCluster(W); clusterExport(cl, ls()); invisible(clusterEvalQ(cl, suppressPackageStartupMessages(library(CompQuadForm))))
rej <- function(v) mean(v < .05, na.rm = TRUE)
# ---- (A)(C)(E)(F) size ----
sz <- do.call(rbind, lapply(c(CAL, HELD), function(d) { o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one_size(r, d), error = function(e) NULL)))
  cat(sprintf("size %s n=%d (%.1f min)\n", d$name, nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))
write.csv(sz, file.path(OUT, "E12_size_perrep.csv"), row.names = FALSE)
agg <- aggregate(cbind(hl, edge, e2, e3, e2_noinf, tau1, tau2, tau3, tau4) ~ cell, sz, rej)
extra <- aggregate(cbind(kappa, rho_hat, deg_final, est_sd, true_sd, tdeg1, tdeg2, tdeg3, tdeg4) ~ cell, sz, mean)
agg <- merge(agg, extra, by = "cell"); agg$set <- ifelse(agg$cell %in% sapply(CAL, `[[`, "name"), "calibration", "held-out"); agg$R <- R
write.csv(agg, file.path(OUT, "E12_size_summary.csv"), row.names = FALSE)
options(width = 250); cat("\n=== (A) size under the final rule + (C) tau sensitivity + (E) scale recovery + (F) uninflated EDGE-2 ===\n")
print(agg[order(agg$set, agg$cell), c("set", "cell", "kappa", "rho_hat", "deg_final", "hl", "edge", "e2", "e3", "e2_noinf", "tau1", "tau2", "tau3", "tau4", "est_sd", "true_sd")], digits = 3, row.names = FALSE)
# ---- (B) power ----
PCELLS <- list(mkA("A_lam100", 100), mkB("B_k10", 400, 40, 416), mkB("B_k25_lam416", 400, 100, 416))
PD <- rbind(data.frame(deg = 2, rel = c(1, 1.5)), data.frame(deg = 3, rel = c(1, 1.5)))
pw <- do.call(rbind, lapply(PCELLS, function(d) do.call(rbind, lapply(seq_len(nrow(PD)), function(k) {
  o <- do.call(rbind, parLapply(cl, seq_len(R), function(r) tryCatch(one_pow(r, d, PD$deg[k], PD$rel[k]), error = function(e) NULL)))
  cat(sprintf("power %s deg%d rel%.1f n=%d (%.1f min)\n", d$name, PD$deg[k], PD$rel[k], nrow(o), as.numeric(difftime(Sys.time(), t0, units = "mins")))); o }))))
write.csv(pw, file.path(OUT, "E12_power_perrep.csv"), row.names = FALSE)
pag <- aggregate(cbind(hl, e2, e3, edge, edge_old) ~ cell + dep_deg + rel, pw, rej)
pag$rho_hat <- round(aggregate(rho_hat ~ cell + dep_deg + rel, pw, mean)$rho_hat, 3)
pag$share3 <- aggregate(deg_final ~ cell + dep_deg + rel, pw, function(v) mean(v == 3))$deg_final
write.csv(pag, file.path(OUT, "E12_power_summary.csv"), row.names = FALSE)
cat("\n=== (B) power by basis: does the kappa clause rescue the cubic at fixed p? (edge = new rule, edge_old = rho2 floor-2 rule) ===\n")
print(pag[order(pag$cell, pag$dep_deg, pag$rel), ], digits = 3, row.names = FALSE)
# ---- (D) oracle Sigma_x ----
ORC <- list(mkB("B_k25_lam416", 400, 100, 416), mkB("B_k40", 400, 160, 416), mkB("H5_k25_ar03", 400, 100, 416, ar = 0.3))
orc <- do.call(rbind, lapply(ORC, function(d) do.call(rbind, parLapply(cl, seq_len(min(R, 200)), function(r) tryCatch(one_oracle(r, d), error = function(e) NULL)))))
stopCluster(cl); write.csv(orc, file.path(OUT, "E12_oracle_perrep.csv"), row.names = FALSE)
oag <- aggregate(cbind(a_hat, a_or, s_hat, s_or, rho_hat, rho_or) ~ cell, orc, mean); write.csv(oag, file.path(OUT, "E12_oracle_summary.csv"), row.names = FALSE)
cat("\n=== (D) sample-covariance plug-in vs oracle Sigma_x ===\n"); print(oag, digits = 3, row.names = FALSE)
cat(sprintf("\ntotal %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
