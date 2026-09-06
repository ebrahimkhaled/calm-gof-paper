# =====================================================================
# figures.R -- the four figures of the CALM paper, built from results/*.csv.
# EK house style (_ek_theme.R): serif, Okabe-Ito, colourblind- and greyscale-safe,
# direct labels, vector PDF at the OUP text width.  Each figure ends in stopifnot()
# guards asserting the claim its caption makes.
# Run from anywhere:  Rscript paper/R/figures.R
# =====================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
HERE <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
ROOT <- dirname(dirname(HERE)); RES <- file.path(ROOT, "results"); FIG <- file.path(ROOT, "paper", "Fig")
source(file.path(HERE, "_ek_theme.R"))
OUP_W <- 446.70827 / 72.27                       # OUP text width in inches; build at true size, never scale down
OI <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73", red = "#D55E00", purple = "#CC79A7", sky = "#56B4E9", grey = "#666666")
cell_lab <- c(A_n500_p5_lam100 = "A, λ=100", A_n500_p5_lam137 = "A, λ=137", A_n500_p5_lam200 = "A, λ=200",
              B_n400_p40_lam416 = "κ=.10", B_n400_p100_lam50 = "κ=.25, λ=50", B_n400_p100_lam416 = "κ=.25, λ=416",
              B_n400_p100_lam1000 = "κ=.25, λ=1000", B_n400_p160_lam416 = "κ=.40")
cell_order <- names(cell_lab)
save_fig <- function(p, name, h) ggsave(file.path(FIG, name), p, width = OUP_W, height = h, device = cairo_pdf)

# ---- Figure 1: the fitted variance is the mechanism ----------------------------------------
e3 <- read.csv(file.path(RES, "E3_jacobian_perrep_R0toR3_run1.csv"))
d1 <- do.call(rbind, lapply(split(e3, e3$cell), function(s) data.frame(cell = s$cell[1],
  R0 = mean(s$S_dec) / mean(s$R0_df_dec), R2 = mean(s$S_dec) / mean(s$R2_df_dec),
  df0 = mean(s$R0_df_dec), df2 = mean(s$R2_df_dec), S = mean(s$S_dec))))
d1$lab <- factor(cell_lab[d1$cell], levels = rev(cell_lab[cell_order]))
pA <- ggplot(d1) + geom_vline(xintercept = 1, linetype = "22", colour = "grey55", linewidth = .4) +
  geom_segment(aes(x = R0, xend = R2, y = lab, yend = lab), colour = "grey70", linewidth = 1.2) +
  geom_point(aes(x = R0, y = lab), colour = OI["grey"], size = 2.6, shape = 21, fill = "white", stroke = 1.1) +
  geom_point(aes(x = R2, y = lab), colour = OI["green"], size = 2.8) +
  annotate("text", x = d1$R0[d1$cell == "A_n500_p5_lam100"] - .012, y = "A, λ=100", label = "first order", hjust = 1, size = 3, colour = OI["grey"], family = "serif") +
  annotate("text", x = d1$R2[d1$cell == "A_n500_p5_lam100"] + .012, y = "A, λ=100", label = "variance at π(β̃)", hjust = 0, size = 3, colour = OI["green"], family = "serif") +
  scale_x_continuous(limits = c(.75, 1.25), breaks = seq(.8, 1.2, .1)) + labs(x = "E[S] / implied E[S]", y = NULL, title = "A") + theme_ek()
b <- d1[d1$cell == "B_n400_p100_lam416", ]
d1b <- data.frame(what = factor(c("first-order
reference", "variance
at π(β̃)", "observed
mean S"), levels = c("first-order
reference", "variance
at π(β̃)", "observed
mean S")), v = c(b$df0, b$df2, b$S))
pB <- ggplot(d1b, aes(x = what, y = v)) + geom_col(fill = unname(c(OI["grey"], OI["green"], "grey30")), width = .6) +
  geom_text(aes(label = sprintf("%.2f", v)), vjust = -.4, size = 3.2, family = "serif") + scale_y_continuous(limits = c(0, 7), expand = c(0, 0)) +
  labs(x = NULL, y = "E[S] under a correct model", title = "B") + theme_ek() + theme(axis.text.x = element_text(size = 8))
save_fig(pA + pB + plot_layout(widths = c(1.35, 1)), "fig1_mechanism.pdf", 3.3)
stopifnot(all(abs(d1$R2[grepl("^A_", d1$cell)] - 1) < .03), all(d1$R0 < .95), b$df0 > 6 & b$df2 < 4.6 & b$S > 4.7 & b$S < 5.0)

# ---- Figure 2: the selection non-centrality is a cubic and is O(kappa) ---------------------
e3c <- read.csv(file.path(RES, "E3c_grouping_mean_perrep.csv"))
vr <- colMeans(e3c[, grep("^v_rand", names(e3c))]); vf <- colMeans(e3c[, grep("^v_fix", names(e3c))])
d2a <- rbind(data.frame(g = 1:10, m = vr, grouping = "on the fitted index (the test)"), data.frame(g = 1:10, m = vf, grouping = "fixed, a function of X only"))
pA2 <- ggplot(d2a, aes(g, m, colour = grouping, shape = grouping)) + geom_hline(yintercept = 0, colour = "grey55", linewidth = .4) +
  geom_line(linewidth = .7) + geom_point(size = 2.4) + scale_colour_manual(values = unname(c(OI["red"], OI["blue"])), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
  scale_x_continuous(breaks = 1:10) + labs(x = "decile of the fitted index", y = "mean corrected residual  E[v_g]", title = "A") +
  theme_ek() + theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8), legend.margin = margin(0, 0, 0, 0))
dl <- read.csv(file.path(RES, "delta_law_cells.csv")); de <- read.csv(file.path(RES, "delta_law_edge_cells.csv"))
d2b <- rbind(data.frame(kappa = dl$kappa, g = dl$meanv2_tr, basis = "decile (all 10 directions)"), data.frame(kappa = de$kappa, g = de$g_edge3, basis = "EDGE-3 (the cubic subspace)"))
pB2 <- ggplot(d2b, aes(kappa, g, colour = basis, shape = basis)) + geom_abline(slope = .24, intercept = 0, colour = OI["grey"], linetype = "22", linewidth = .5) +
  geom_abline(slope = 1.15, intercept = 0, colour = OI["red"], linetype = "22", linewidth = .5) + geom_point(size = 2.6) +
  scale_colour_manual(values = unname(c(OI["grey"], OI["red"])), name = NULL) + scale_shape_manual(values = c(16, 15), name = NULL) +
  annotate("text", x = .33, y = .24 * .33 + .03, label = "0.24 κ", size = 3, family = "serif", colour = OI["grey"]) +
  annotate("text", x = .30, y = 1.15 * .30 + .05, label = "1.15 κ", size = 3, family = "serif", colour = OI["red"]) +
  scale_x_continuous(limits = c(0, .42), breaks = c(0, .1, .25, .4)) + labs(x = "κ = p / n", y = "||E v||² / tr(Σ)", title = "B") +
  theme_ek() + theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 8), legend.margin = margin(0, 0, 0, 0))
save_fig(pA2 + pB2, "fig2_selection.pdf", 3.6)
stopifnot(sum(vr^2) > 3 * sum(vf^2), cor(dl$kappa, dl$meanv2_tr) > .9, cor(de$kappa, de$g_edge3) > .9, all(de$share_in_Z3[de$kappa > .05] > .85))

# ---- Figure 3: the attenuation law as a design rule --------------------------------------------
e8 <- read.csv(file.path(RES, "E8_basis_degree_summary.csv"))
e8 <- e8[e8$deg > 0 & e8$rel == 1.0, ]
long <- do.call(rbind, lapply(c("p_dec", "p_e1", "p_e2", "p_e3"), function(v) data.frame(cell = e8$cell, deg = e8$deg, rho = e8$rho_hat, basis = v, power = e8[[v]])))
long$basis <- factor(c(p_dec = "deciles", p_e1 = "EDGE-1", p_e2 = "EDGE-2", p_e3 = "EDGE-3")[long$basis], levels = c("deciles", "EDGE-1", "EDGE-2", "EDGE-3"))
long$cell <- factor(c(A_p5 = "κ = 0.01  (ρ̂ ≈ .94)", B_p40 = "κ = 0.10  (ρ̂ ≈ .71)", B_p100 = "κ = 0.25  (ρ̂ ≈ .69)")[long$cell], levels = c("κ = 0.01  (ρ̂ ≈ .94)", "κ = 0.10  (ρ̂ ≈ .71)", "κ = 0.25  (ρ̂ ≈ .69)"))
long$departure <- factor(ifelse(long$deg == 2, "degree-2 departure", "degree-3 departure"))
p3 <- ggplot(long, aes(basis, power, fill = departure)) + geom_col(position = position_dodge(width = .7), width = .65) + ek_nominal(.05) +
  facet_wrap(~cell, nrow = 1) + scale_fill_manual(values = unname(c(OI["orange"], OI["blue"])), name = NULL) + scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "power at the 5% level", title = NULL) +
  theme_ek() + theme(legend.position = "bottom", axis.text.x = element_text(size = 8))
save_fig(p3, "fig3_attenuation.pdf", 3.0)
w <- function(cell, deg, v) e8[[v]][e8$cell == cell & e8$deg == deg]
stopifnot(w("B_p100", 2, "p_e2") > w("B_p100", 2, "p_e3"), w("B_p100", 2, "p_e2") > w("B_p100", 2, "p_dec"), w("A_p5", 3, "p_e3") > .5, w("B_p100", 3, "p_e3") < .10, w("B_p40", 3, "p_e3") < .10)

# ---- Figure 4: size of the references across the eight cells -----------------------------------
e2b <- read.csv(file.path(RES, "E2b_deflation_summary.csv")); e2b <- e2b[e2b$grouping == "random", ]
e3b <- read.csv(file.path(RES, "E3b_denoised_summary.csv")); e9 <- read.csv(file.path(RES, "E9_kappa_inflation_summary.csv")); e9b <- read.csv(file.path(RES, "E9b_edge_inflation_summary.csv"))
mk <- function(cell, ref, basis, size) data.frame(cell = cell, ref = ref, basis = basis, size = size)
d4 <- rbind(mk(e2b$cell, "first order", "SC.HL", e2b$rej_dec), mk(e2b$cell, "first order", "SC.EDGE", e2b$rej_edge),
            mk(e3b$cell[e3b$ref == "R2"], "variance at π(β̃)", "SC.HL", e3b$size_dec[e3b$ref == "R2"]), mk(e3b$cell[e3b$ref == "R2"], "variance at π(β̃)", "SC.EDGE", e3b$size_edge[e3b$ref == "R2"]),
            mk(e3b$cell[e3b$ref == "R5"], "de-noised variance", "SC.HL", e3b$size_dec[e3b$ref == "R5"]), mk(e3b$cell[e3b$ref == "R5"], "de-noised variance", "SC.EDGE", e3b$size_edge[e3b$ref == "R5"]),
            { e4c3 <- read.csv(file.path(RES, "E4c3_size_summary.csv"))   # FINAL default (degree rule rho2, 2026-09-04): the source of Table 1's CALM columns
              e4c3$cell <- sub("^A_p", "A_n500_p", sub("^B_p", "B_n400_p", e4c3$cell))
              rbind(mk(e4c3$cell, "CALM", "SC.HL", e4c3$hl), mk(e4c3$cell, "CALM", "SC.EDGE", e4c3$edge_k)) })
d4$ref <- factor(d4$ref, levels = c("first order", "variance at π(β̃)", "de-noised variance", "CALM")); d4$lab <- factor(cell_lab[d4$cell], levels = cell_lab[cell_order])
d4$basis <- factor(d4$basis, levels = c("SC.HL", "SC.EDGE"))
se <- sqrt(.05 * .95 / 500)
p4 <- ggplot(d4, aes(lab, size, colour = ref, shape = ref)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = .05 - 2 * se, ymax = .05 + 2 * se, fill = "grey88", alpha = .6) + ek_nominal(.05) +
  geom_point(position = position_dodge(width = .6), size = 2.4) + facet_wrap(~basis, ncol = 1) +
  scale_colour_manual(values = unname(c(OI["grey"], OI["sky"], OI["orange"], OI["green"])), name = NULL) + scale_shape_manual(values = c(21, 17, 15, 16), name = NULL) +
  scale_y_continuous(limits = c(0, .14), breaks = c(0, .05, .10)) + labs(x = NULL, y = "rejection rate of a correct model (nominal .05)", title = NULL) +
  theme_ek() + theme(legend.position = "bottom", axis.text.x = element_text(size = 8, angle = 20, hjust = 1))
save_fig(p4, "fig4_size.pdf", 4.6)
calm <- d4[d4$ref == "CALM", ]; stopifnot(all(calm$size <= .065), all(d4$size[d4$ref == "first order" & grepl("^B_", d4$cell)] < .04))
cat("figures written:", paste(list.files(FIG, pattern = "pdf$"), collapse = ", "), "\n")
