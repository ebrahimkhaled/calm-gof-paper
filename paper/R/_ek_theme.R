## _ek_theme.R — one shared, colourblind-safe (Okabe-Ito), greyscale-legible figure system
## for the EDGE paper. source() at the top of every plotting script. Display keys are the
## paper's names (EDGE-poly2/3, EDGE-stk); plotting scripts recode CSV test names
## (DEF.poly2 -> EDGE-poly2, HL-equalwidth -> HLw, Pigeon-Heyse -> PH, DEF.stukel -> EDGE-stk).
## base_family="sans" => Arial on Windows (Springer Helvetica/Arial requirement, A8).
suppressMessages({ library(ggplot2); library(scales) })

ek_pal <- c("EDGE-poly2" = "#E69F00", "EDGE-poly3" = "#D55E00", "EDGE-stk" = "#CC79A7",
            "EF" = "#0072B2", "HL" = "#56B4E9", "HLw" = "#80B1D3", "PH" = "#009E73",
            "Tsiatis" = "#7F7F7F", "Stukel" = "#000000", "Xie" = "#B4A0C7", "PR" = "#BBBBBB")
ek_lty <- c("EDGE-poly2" = "solid", "EDGE-poly3" = "solid", "EDGE-stk" = "solid", "Stukel" = "solid",
            "EF" = "22", "HL" = "42", "HLw" = "42", "PH" = "42", "Tsiatis" = "42", "Xie" = "42", "PR" = "42")
ek_levels <- c("EDGE-poly3", "EDGE-poly2", "EDGE-stk", "Stukel", "EF", "HL", "HLw", "PH", "Tsiatis", "Xie", "PR")

## map raw CSV Test names -> display keys above
ek_relabel <- function(x) {
  m <- c("DEF.poly2" = "EDGE-poly2", "DEF.poly3" = "EDGE-poly3", "DEF.stukel" = "EDGE-stk",
         "HL-equalwidth" = "HLw", "Pigeon-Heyse" = "PH", "Pulkstenis-Robinson" = "PR")
  x <- as.character(x); ifelse(x %in% names(m), m[x], x)
}
ek_factor <- function(x) factor(ek_relabel(x), levels = ek_levels)

scale_color_ek <- function(...) scale_color_manual(values = ek_pal, breaks = ek_levels, name = NULL, ...)
scale_lty_ek   <- function(...) scale_linetype_manual(values = ek_lty, breaks = ek_levels, guide = "none", ...)
scale_fill_alpha_ek <- function(mid = 0.05, lim = c(0, 0.12))
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = mid,
                       limits = lim, oob = scales::squish, name = expression(hat(alpha)))

theme_ek <- function(base_size = 10.5, base_family = "serif") theme_minimal(base_size, base_family) %+replace%
  theme(plot.title = element_text(face = "plain", size = rel(1.05), hjust = 0.5, margin = margin(b = 4)),
        plot.subtitle = element_text(size = rel(0.82), hjust = 0.5, colour = "grey35"),
        plot.title.position = "plot", panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "plain", size = rel(0.9)), legend.position = "bottom")

ek_nominal <- function(a = 0.05) geom_hline(yintercept = a, linetype = "13", colour = "grey55", linewidth = 0.4)
ek_diag    <- function() geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey55", linewidth = 0.4)

## Springer figure widths (mm): single-column text area 174, double-column 84 (A8).
EK_W1 <- 84 / 25.4; EK_W2 <- 174 / 25.4; EK_DPI <- 320
