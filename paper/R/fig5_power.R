# =====================================================================
# fig5_power.R -- Figure 5 (power of CALM against the bootstrap and the rivals) + results/power_master_summary.csv
# Sources, all on IDENTICAL datasets (SEED + rep):  E4c2 (CALM variants; the FINAL design is re-derived here with the floor-2 rule),
#   E4b (null) + E4c (departures) for prepivoting, E4d for the rivals (GRP, classical on the MLE, CV-HL, calibration slope, Spiegelhalter,
#   USS, uncorrected HL).  Tests that do not hold the level are shown greyed in the figure and get size-adjusted power in the table.
# Run from anywhere:  Rscript paper/R/fig5_power.R
# =====================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(dirname(HERE)); RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "paper", "Fig")
source(file.path(HERE, "_ek_theme.R"))
OUP_W <- 446.70827 / 72.27
OI <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73", red = "#D55E00", purple = "#CC79A7", sky = "#56B4E9", grey = "#666666")
save_fig <- function(p, name, h) ggsave(file.path(FIG, name), p, width = OUP_W, height = h, device = cairo_pdf)
thr <- 0.2^(1 / 6)
key <- c("rep", "cell", "dep", "rel")

# ---- CALM (final design) from E4c2 per rep ----
v <- read.csv(file.path(RES, "E4c2_variants_perrep.csv"))
calm <- data.frame(v[key], calm_hl = v$hl_R5_inf, calm_edge = ifelse(v$rho_hat >= thr, v$e3_R5_inf, v$e2_R5_inf), rho_hat = v$rho_hat)
# ---- prepivoting: E4b (null rows) + E4c (departures) ----
e4b <- read.csv(file.path(RES, "E4b_final_perrep.csv")); e4b <- e4b[e4b$dep == "null", c(key, "pp_hl", "pp_edge")]   # E4b's own departures (cubic, small x1x2) belong to the boundary section
pp <- rbind(e4b, do.call(rbind, lapply(list.files(RES, "^E4c_perrep_", full.names = TRUE), function(f) read.csv(f)[, c(key, "pp_hl", "pp_edge")])))
# ---- rivals: E4d ----
rv <- if (file.exists(file.path(RES, "E4d_rivals_perrep.csv"))) read.csv(file.path(RES, "E4d_rivals_perrep.csv")) else
      do.call(rbind, lapply(list.files(RES, "^E4d_perrep_", full.names = TRUE), read.csv))   # partial run: per-cell files
m <- merge(merge(calm, pp, by = key, all = TRUE), rv, by = key, all = TRUE)
m$cell <- factor(m$cell, levels = c("A_p5", "B_p40", "B_p100", "B_p160"), labels = c("fixed p (κ = 0.01)", "κ = 0.10", "κ = 0.25", "κ = 0.40"))
tests <- c(calm_hl = "CALM, SC.HL", calm_edge = "CALM, SC.EDGE", pp_hl = "prepivot, SC.HL", pp_edge = "prepivot, SC.EDGE", grp = "GRP (Janková et al.)",
           mle_hl = "HL on the MLE fit", mle_edge = "EDGE on the MLE fit", cvhl = "HL on CV predictions",
           calslope = "calibration slope test", spieg = "Spiegelhalter z", uss = "unweighted sum of squares", unc_hl = "uncorrected HL")
rej <- function(p) mean(p < .05, na.rm = TRUE)
agg <- aggregate(m[, names(tests)], by = m[, c("cell", "dep", "rel")], FUN = rej)
agg <- merge(agg, aggregate(mle_ok ~ cell + dep + rel, m, function(z) mean(z, na.rm = TRUE)), by = c("cell", "dep", "rel"), all.x = TRUE)   # share of datasets where the MLE existed
# size-adjusted power: critical value = empirical 5% quantile of the same test's null p-values in the same cell
adj <- do.call(rbind, lapply(split(m, m$cell), function(dc) { nul <- dc[dc$dep == "null", ]
  do.call(rbind, lapply(split(dc, list(dc$dep, dc$rel), drop = TRUE), function(dk) { o <- data.frame(cell = dk$cell[1], dep = dk$dep[1], rel = dk$rel[1])
    for (t in names(tests)) { cv <- suppressWarnings(quantile(nul[[t]], .05, na.rm = TRUE))
      # when the null p-values have underflown (calibration slope, Spiegelhalter, USS under strong shrinkage) the ranking is lost and
      # size-adjusted power is not computable from p-values: report NA rather than a degenerate 0 or 1
      o[[paste0(t, "_adj")]] <- if (is.finite(cv) && cv > 1e-10) mean(dk[[t]] <= cv, na.rm = TRUE) else NA }; o })) }))
agg <- merge(agg, adj, by = c("cell", "dep", "rel"))
write.csv(agg, file.path(RES, "power_master_summary.csv"), row.names = FALSE)
# which tests hold the level: max null rejection over the four cells <= 0.075 (2.5 MC s.e. at R = 300-500)
nul <- agg[agg$dep == "null", names(tests)]; valid <- names(tests)[sapply(nul, function(z) all(z <= .075, na.rm = TRUE))]
cat("tests holding the level in every cell:", paste(tests[valid], collapse = " | "), "\n")
cat("null rejection rates by cell (all tests):\n"); print(cbind(cell = as.character(agg$cell[agg$dep == "null"]), round(nul, 3)), row.names = FALSE)

# ---- Figure 5: power curves, valid tests solid, invalid tests greyed/dashed (raw rejection rate) ----
long <- do.call(rbind, lapply(names(tests), function(t) data.frame(cell = agg$cell, dep = agg$dep, rel = agg$rel, test = tests[[t]], power = agg[[t]], valid = t %in% valid)))
long <- long[!is.na(long$power), ]
# per-cell validity: a test's power curve is drawn in a cell only where its size in that cell is <= 0.075; elsewhere it goes to the table
# (size-adjusted power) -- a rejection rate above a broken level is not power.
sz <- long[long$dep == "null", c("cell", "test", "power")]; names(sz)[3] <- "size_cell"
long <- merge(long, sz, by = c("cell", "test")); long <- long[!is.na(long$size_cell) & long$size_cell <= .075, ]
add0 <- function(d, depname) { z <- d[d$dep == "null", ]; z$dep <- depname; rbind(d[d$dep == depname, ], z) }
lq <- add0(long, "quadratic"); lo <- add0(long, "offindex")
show <- c("CALM, SC.HL", "CALM, SC.EDGE", "prepivot, SC.HL", "prepivot, SC.EDGE", "GRP (Janková et al.)", "HL on the MLE fit")   # CV-HL never holds its level -> table only
cols <- c("CALM, SC.HL" = OI[["green"]], "CALM, SC.EDGE" = OI[["blue"]], "prepivot, SC.HL" = OI[["orange"]], "prepivot, SC.EDGE" = OI[["red"]],
          "GRP (Janková et al.)" = OI[["purple"]], "HL on the MLE fit" = "grey35", "HL on CV predictions" = OI[["sky"]])
shp <- c("CALM, SC.HL" = 16, "CALM, SC.EDGE" = 15, "prepivot, SC.HL" = 1, "prepivot, SC.EDGE" = 0, "GRP (Janková et al.)" = 17, "HL on the MLE fit" = 4, "HL on CV predictions" = 3)
pl <- function(d, ttl, xl, legend) { d <- d[d$test %in% show, ]; d$test <- factor(d$test, levels = show)
  ggplot(d, aes(rel, power, colour = test, shape = test, linetype = valid)) + ek_nominal(.05) + geom_line(linewidth = .55) + geom_point(size = 2) +
    facet_wrap(~cell, nrow = 1) + scale_colour_manual(values = cols, name = NULL, drop = FALSE) + scale_shape_manual(values = shp, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), guide = "none") + scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
    labs(x = xl, y = "rejection rate at the 5% level", title = ttl) + theme_ek() +
    theme(legend.position = if (legend) "bottom" else "none", legend.text = element_text(size = 8)) + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) }
p5 <- pl(lq, "A. Departure along the index: Hermite quadratic", NULL, FALSE) /
      pl(lo, "B. Departure off the index: omitted x1·x2", "relative size of the departure (0 = size)", TRUE)
save_fig(p5, "fig5_power.pdf", 6.2)
cat("figure written: fig5_power.pdf\n")
# ---- LaTeX rows for the power table (valid tests, raw rates) ----
fmt <- function(x) ifelse(is.na(x), "---", sprintf("%.3f", x))
tab <- agg[order(agg$cell, agg$dep, agg$rel), ]
rows <- apply(tab, 1, function(r) sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", r[["cell"]], r[["dep"]], r[["rel"]],
        fmt(as.numeric(r[["calm_hl"]])), fmt(as.numeric(r[["calm_edge"]])), fmt(as.numeric(r[["pp_hl"]])), fmt(as.numeric(r[["pp_edge"]])), fmt(as.numeric(r[["grp"]])), fmt(as.numeric(r[["mle_hl"]])), fmt(as.numeric(r[["cvhl"]]))))
writeLines(rows, file.path(RES, "power_table_rows.tex")); cat(rows, sep = "\n")
