# Evaluate the FLOOR-2 hybrid degree rule offline: k = 3 if rho_hat^6 >= 0.2 (rho_hat >= 0.765), else k = 2 (never 1).
# Uses E4c2 (variants per rep, has rho_hat, e2_R5_inf, e3_R5_inf) + prepivot EDGE from E4b (null) and E4c (departures); E4c3 for the 8 size cells.
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1]))), "..", "results"))
thr <- 0.2^(1 / 6)
d <- read.csv("E4c2_variants_perrep.csv"); d$hyb <- ifelse(d$rho_hat >= thr, d$e3_R5_inf, d$e2_R5_inf)
pp <- rbind(read.csv("E4b_final_perrep.csv")[, c("rep", "cell", "dep", "rel", "pp_edge")],
            do.call(rbind, lapply(list.files(".", "^E4c_perrep_", full.names = TRUE), function(f) read.csv(f)[, c("rep", "cell", "dep", "rel", "pp_edge")])))
m <- merge(d, pp, by = c("rep", "cell", "dep", "rel"))
r <- function(v) mean(v < .05, na.rm = TRUE)
a <- aggregate(cbind(hyb, e2_R5_inf, e3_R5_inf, ad_R5_inf, pp_edge) ~ cell + dep + rel, m, r)
a$share_k3 <- aggregate(I(rho_hat >= thr) ~ cell + dep + rel, m, mean)[[4]]
options(width = 200); cat(sprintf("threshold rho_hat >= %.3f\n", thr)); print(a[order(a$cell, a$dep, a$rel), ], digits = 3, row.names = FALSE)
s <- read.csv("E4c3_size_perrep.csv"); s$hyb <- ifelse(s$rho_hat >= thr, s$e3, s$edge_k)
cat("\nE4c3 sizes (8 cells):\n"); print(aggregate(cbind(hl, hyb, edge_k, e3) ~ cell, s, r), digits = 3, row.names = FALSE)
write.csv(a, "hybrid_rule_eval.csv", row.names = FALSE)
