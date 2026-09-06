# EDGE-specific non-centrality fraction: ||P_Z E[v]||^2 / df_edge, Z ~ poly(group index, 3) (equal-frequency groups)
e3 <- read.csv("results/E3_jacobian_perrep_R0toR3_run1.csv"); e3b <- read.csv("results/E3b_denoised_perrep.csv"); e9 <- read.csv("results/E9_kappa_inflation_perrep.csv")
Z3 <- as.matrix(poly(1:10, 3)); Z2 <- as.matrix(poly(1:10, 2)); P <- function(Z) Z %*% solve(crossprod(Z), t(Z))
rows <- list()
for (cn in unique(e3$cell)) { s <- e3[e3$cell == cn, ]; vbar <- colMeans(s[, grep("^v[0-9]+$", names(s))]); b <- e3b[e3b$cell == cn, ]; k <- e9[e9$cell == cn, ]
  parts <- strsplit(cn, "_")[[1]]; n <- as.integer(sub("n", "", parts[2])); p <- as.integer(sub("p", "", parts[3]))
  rows[[length(rows) + 1]] <- data.frame(cell = cn, kappa = p / n, hbar = mean(k$hbar), rho = mean(k$rho),
    g_dec = sum(vbar^2) / mean(b$R5_df_dec), g_edge3 = sum((P(Z3) %*% vbar)^2) / mean(b$R5_df_edge),
    g_edge2 = sum((P(Z2) %*% vbar)^2) / (mean(b$R5_df_edge) * 2 / 3), share_in_Z3 = sum((P(Z3) %*% vbar)^2) / sum(vbar^2)) }
d <- do.call(rbind, rows); options(width = 200); print(d[order(d$kappa), ], row.names = FALSE, digits = 3)
cat("\nlaws through origin for the EDGE-3 fraction:\n")
for (f in c("kappa", "I(sqrt(hbar))", "I(1-rho^2)")) { m <- lm(as.formula(paste("g_edge3 ~ 0 +", f)), d); cat(sprintf("  %-14s coef %.3f  R2 %.3f  max|res| %.3f\n", f, coef(m), 1 - sum(resid(m)^2) / sum(d$g_edge3^2), max(abs(resid(m))))) }
m <- lm(g_edge2 ~ 0 + kappa, d); cat(sprintf("  EDGE-2 fraction ~ kappa: coef %.3f  R2 %.3f\n", coef(m), 1 - sum(resid(m)^2) / sum(d$g_edge2^2)))
write.csv(d, "results/delta_law_edge_cells.csv", row.names = FALSE)
