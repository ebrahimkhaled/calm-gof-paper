args <- commandArgs(trailingOnly = TRUE); f <- if (length(args)) args[1] else "results/E3_jacobian_perrep.csv"
res <- read.csv(f); options(width = 200)
refs <- intersect(c("R0","R1","R2","R3","R4"), sub("_df_dec$", "", grep("_df_dec$", names(res), value = TRUE)))
cat("rows:", nrow(res), " cells:", length(unique(res$cell)), " refs:", paste(refs, collapse = ","), "\n")
rows <- list()
for (cn in unique(res$cell)) {
  s <- res[res$cell == cn, ]; vbar <- colMeans(s[, grep("^v[0-9]+$", names(s)), drop = FALSE])
  for (k in refs) rows[[length(rows) + 1]] <- data.frame(cell = cn, ref = k,
    defl_dec = mean(s$S_dec) / mean(s[[paste0(k, "_df_dec")]]), size_dec = mean(s[[paste0(k, "_p_dec")]] < 0.05),
    defl_edge = mean(s$S_edge, na.rm = TRUE) / mean(s[[paste0(k, "_df_edge")]], na.rm = TRUE), size_edge = mean(s[[paste0(k, "_p_edge")]] < 0.05, na.rm = TRUE),
    implied_dec = mean(s[[paste0(k, "_df_dec")]]), mean_S = mean(s$S_dec), w_ratio = mean(s$w_ratio),
    w_ratio_recal = if ("w_ratio_recal" %in% names(s)) mean(s$w_ratio_recal) else NA, meanv2_tr = sum(vbar^2) / mean(s$tr_OmMLE))
}
a <- do.call(rbind, rows)
for (v in c("defl_dec", "size_dec", "defl_edge", "size_edge", "w_ratio", "w_ratio_recal", "meanv2_tr")) a[[v]] <- sprintf("%.3f", a[[v]])
for (v in c("implied_dec", "mean_S")) a[[v]] <- sprintf("%.2f", a[[v]])
print(a, row.names = FALSE)
