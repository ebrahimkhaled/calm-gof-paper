suppressPackageStartupMessages({library(glmnet); library(TH.data)})
source("R/closedform_engine.R")
data("GlaucomaM", package = "TH.data")
y <- as.numeric(GlaucomaM$Class == "glaucoma")
X <- scale(as.matrix(GlaucomaM[, setdiff(names(GlaucomaM), "Class")])); n <- length(y)
set.seed(1); cv <- cv.glmnet(X, y, family = "binomial", alpha = 0); l1se <- cv$lambda.1se * n
r <- closedform.gof(X, y, lambda = l1se, G = 10, uncorrected = TRUE); print(r)
cat("prepivoted (shrink.gof, CRAN):  SC.HL 0.026   SC.EDGE 0.038   -- Table 5 of the Series C paper: 0.026 / 0.034\n")
cat("\nobservables:\n"); print(round(r$obs, 3))
