source("R/closedform_engine.R"); suppressPackageStartupMessages({library(TH.data); library(ebrahim.gof)})
data("GlaucomaM", package = "TH.data"); y <- as.numeric(GlaucomaM$Class == "glaucoma"); X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")]))
t_r2 <- system.time(for (i in 1:5) closedform.gof(X, y, 295.3, reference = "R2"))[3] / 5
t_r5 <- system.time(for (i in 1:5) closedform.gof(X, y, 295.3, reference = "R5"))[3] / 5
t_pp <- system.time(shrink.gof(X, y, lambda = 295.3, G = 10, B = 499, seed = 1))[3]
cat(sprintf("TIMING glaucoma (n=196, p=62): closed form R2 %.2f s | R5 %.2f s | prepivot B=499 %.1f s  ->  speedup R2 %.0fx, R5 %.0fx\n", t_r2, t_r5, t_pp, t_pp / t_r2, t_pp / t_r5))
