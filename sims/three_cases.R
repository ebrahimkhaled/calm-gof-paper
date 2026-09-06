# =====================================================================
# three_cases.R -- Correct form / wrong form ALONG the index / wrong form
#                  OFF the index, each fitted by ridge, with and without
#                  recalibration.
#
# The question this answers: what does recalibration repair, and what does
# it hide?  Calibration (ECE, slope, held-out HL) and discrimination (AUC)
# can both look fine while the individual probabilities are wrong.  The one
# quantity that tells the truth is the error against the TRUE risk, which
# only a simulation can see.  The structural test (shrink.gof, the paper's
# SC.HL / SC.EDGE) is run on the training data as the H0^str verdict.
#
# Designs are the registered ones of the Series C paper.
#   A: n = 500, p = 5,   independent N(0,1),      lambda = 100  (lambda_g 0.20)
#   B: n = 400, p = 100, AR(0.7) correlation,     lambda = 416  (lambda_g 1.04)
#
# Usage:  Rscript three_cases.R A [R] [B_boot] [workers]
# Output: results/three_cases_<D>_perrep.csv   one row per rep x case x recal
#         results/three_cases_<D>_summary.csv  the table
# =====================================================================
suppressPackageStartupMessages({
  library(ebrahim.gof); library(splines); library(parallel)
})
args    <- commandArgs(trailingOnly = TRUE)
DESIGN  <- if (length(args) >= 1) args[1] else "A"
R       <- if (length(args) >= 2) as.integer(args[2]) else 200
B_BOOT  <- if (length(args) >= 3) as.integer(args[3]) else 199
WORKERS <- if (length(args) >= 4) as.integer(args[4]) else 20
SEED    <- 20260904
stopifnot(DESIGN %in% c("A", "B"))
OUT <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))), "results")
if (!dir.exists(OUT)) dir.create(OUT)

# ---- registered designs -------------------------------------------------
if (DESIGN == "A") {
  n_train <- 500; p <- 5
  beta0   <- c(0.5, -0.4, 0.3, -0.3, 0.2)
  lambda  <- 100
  Sigma_chol <- diag(p)
} else {
  n_train <- 400; p <- 100
  beta0   <- rep(c(0.35, -0.30, 0.25, -0.20, 0.15), 20) * 0.8887
  lambda  <- 416
  Sigma   <- 0.7^abs(outer(1:p, 1:p, "-"))
  Sigma_chol <- chol(Sigma)
}
n_recal <- n_train        # a validation sample of the same size, for recalibration
n_hold  <- 1000           # a second held-out sample, realistic size, for the calibration test
n_test  <- 20000          # a large sample so evaluation noise is negligible

# ---- departures: ONE dial for both, in both designs ----------------------
# Both departures inject the same amount of extra variation into the logit:
# sd(departure) = REL * sd(eta).  REL is the paper's "rel" (Series C, Sec 4.1);
# the index departure stays MONOTONE in eta iff REL < sqrt(6)/3 = 0.8165, so
# "same ordering as the truth" is exactly true, not approximately.
REL <- 0.70
stopifnot(REL < sqrt(6) / 3)
sd_eta  <- if (DESIGN == "A") sqrt(sum(beta0^2)) else sqrt(drop(t(beta0) %*% (0.7^abs(outer(1:p, 1:p, "-"))) %*% beta0))
sd_x1x2 <- if (DESIGN == "A") 1 else sqrt(1 + 0.7^2)   # Var(x1 x2) = 1 + rho^2 for bivariate normal
GAM_OFF <- REL * sd_eta / sd_x1x2                          # so sd(GAM_OFF * x1 * x2) = REL * sd(eta)
GAM_INDEX <- REL                                           # recorded for the table

gen_X <- function(n) matrix(rnorm(n * p), n, p) %*% Sigma_chol
true_logit <- function(X, case) {
  eta <- drop(X %*% beta0)
  switch(case,
    correct  = eta,
    index    = { Z <- eta / sd_eta; eta + sd_eta * REL * (Z^3 - 3 * Z) / sqrt(6) },  # orthogonal (Hermite) cubic: 100% nonlinear, monotone
    offindex = eta + GAM_OFF   * X[, 1] * X[, 2]) # NOT a function of eta
}
expit <- function(z) 1 / (1 + exp(-z))

# ---- ridge logistic by penalized IRLS (matches glmnet to 1e-9, intercept free)
ridge_fit <- function(X, y, lambda, tol = 1e-9, maxit = 60) {
  X1 <- cbind(1, X); D <- diag(c(0, rep(1, ncol(X))))
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta); pi <- expit(eta); w <- pmax(pi * (1 - pi), 1e-8)
    z  <- eta + (y - pi) / w
    bn <- drop(solve(crossprod(X1, w * X1) + lambda * D, crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  beta
}
lp <- function(beta, X) drop(cbind(1, X) %*% beta)

# ---- metrics ---------------------------------------------------------
auc <- function(pred, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(rank(pred)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
ece10 <- function(pred, y) {
  g <- cut(rank(pred, ties.method = "first"), 10, labels = FALSE)
  w <- tabulate(g, 10) / length(y)
  sum(w * abs(tapply(y, g, mean) - tapply(pred, g, mean)))
}
cal_slope_int <- function(pred, y) {
  lg <- qlogis(pmin(pmax(pred, 1e-6), 1 - 1e-6))
  co <- coef(glm(y ~ lg, family = binomial))
  c(int = unname(co[1]), slope = unname(co[2]))
}
hl_heldout <- function(pred, y, G = 10) {
  # model is FIXED given the training data, so on held-out data the
  # reference is chi-square with G df (no estimation effect).
  g  <- cut(rank(pred, ties.method = "first"), G, labels = FALSE)
  O  <- tapply(y, g, sum); E <- tapply(pred, g, sum); V <- tapply(pred * (1 - pred), g, sum)
  pchisq(sum((O - E)^2 / V), df = G, lower.tail = FALSE)
}

# ---- recalibration maps, fitted on the recalibration sample ------------
recal_none     <- function(lp_r, y_r) function(lp_new) expit(lp_new)
recal_linear   <- function(lp_r, y_r) {
  f <- glm(y_r ~ lp_r, family = binomial)
  function(lp_new) predict(f, newdata = data.frame(lp_r = lp_new), type = "response")
}
recal_flexible <- function(lp_r, y_r) {
  f <- glm(y_r ~ ns(lp_r, df = 4), family = binomial)
  function(lp_new) predict(f, newdata = data.frame(lp_r = lp_new), type = "response")
}
RECALS <- list(none = recal_none, linear = recal_linear, flexible = recal_flexible)
CASES  <- c("correct", "index", "offindex")

# ---- one replication ---------------------------------------------------
one_rep <- function(r) {
  set.seed(SEED + r)
  rows <- list()
  for (case in CASES) {
    Xtr <- gen_X(n_train); ptr <- expit(true_logit(Xtr, case)); ytr <- rbinom(n_train, 1, ptr)
    Xre <- gen_X(n_recal); pre <- expit(true_logit(Xre, case)); yre <- rbinom(n_recal, 1, pre)
    Xho <- gen_X(n_hold);  pho <- expit(true_logit(Xho, case)); yho <- rbinom(n_hold,  1, pho)
    Xte <- gen_X(n_test);  pte <- expit(true_logit(Xte, case)); yte <- rbinom(n_test,  1, pte)

    b   <- ridge_fit(Xtr, ytr, lambda)
    lpr <- lp(b, Xre); lph <- lp(b, Xho); lpt <- lp(b, Xte)

    # the structural test, on the TRAINING data (H0^str)
    sg <- tryCatch(shrink.gof(Xtr, ytr, lambda = lambda, G = 10,
                              basis = c("edge", "decile"), B = B_BOOT, seed = SEED + r),
                   error = function(e) NULL)
    p_edge <- if (is.null(sg)) NA else sg$SC.EDGE$p.value
    p_hl   <- if (is.null(sg)) NA else sg$SC.HL$p.value

    # oracle row: the true probabilities themselves
    cs <- cal_slope_int(pte, yte)
    rows[[length(rows) + 1]] <- data.frame(rep = r, case = case, recal = "oracle",
      ece = ece10(pte, yte), slope = cs["slope"], intercept = cs["int"],
      brier = mean((pte - yte)^2), mae_true = 0, auc = auc(pte, yte),
      hl_heldout_p = hl_heldout(pho, yho), sc_edge_p = NA, sc_hl_p = NA)

    for (rn in names(RECALS)) {
      map <- RECALS[[rn]](lpr, yre)
      pt  <- pmin(pmax(map(lpt), 1e-6), 1 - 1e-6)
      ph  <- pmin(pmax(map(lph), 1e-6), 1 - 1e-6)
      cs  <- cal_slope_int(pt, yte)
      rows[[length(rows) + 1]] <- data.frame(rep = r, case = case, recal = rn,
        ece = ece10(pt, yte), slope = cs["slope"], intercept = cs["int"],
        brier = mean((pt - yte)^2), mae_true = mean(abs(pt - pte)), auc = auc(pt, yte),
        hl_heldout_p = hl_heldout(ph, yho), sc_edge_p = p_edge, sc_hl_p = p_hl)
    }
  }
  do.call(rbind, rows)
}

# ---- run ---------------------------------------------------------------
options(width = 220); t0 <- Sys.time()
cl <- makeCluster(min(WORKERS, R))
clusterExport(cl, setdiff(ls(), c("cl", "t0")))
invisible(clusterEvalQ(cl, suppressPackageStartupMessages({library(ebrahim.gof); library(splines)})))
res <- do.call(rbind, parLapply(cl, seq_len(R), one_rep))
stopCluster(cl)
rownames(res) <- NULL
write.csv(res, file.path(OUT, sprintf("three_cases_%s_perrep.csv", DESIGN)), row.names = FALSE)

# ---- summary table -----------------------------------------------------
res$recal <- factor(res$recal, levels = c("oracle", "none", "linear", "flexible"))
res$case  <- factor(res$case,  levels = CASES)
summ <- aggregate(cbind(ece, slope, brier, mae_true, auc) ~ case + recal, data = res, FUN = mean)
rej  <- aggregate(cbind(hl_heldout_rej = hl_heldout_p < 0.05,
                        sc_edge_rej = sc_edge_p < 0.05,
                        sc_hl_rej   = sc_hl_p   < 0.05) ~ case + recal, data = res,
                  FUN = function(v) mean(v, na.rm = TRUE), na.action = na.pass)
summ <- merge(summ, rej, by = c("case", "recal"))
summ <- summ[order(summ$case, summ$recal), ]
summ$design <- DESIGN; summ$R <- R; summ$B_boot <- B_BOOT
summ$rel <- REL; summ$gam_off <- GAM_OFF; summ$sd_eta <- sd_eta
write.csv(summ, file.path(OUT, sprintf("three_cases_%s_summary.csv", DESIGN)), row.names = FALSE)

cat(sprintf("\n=== DESIGN %s | R = %d | B_boot = %d | %.1f min ===\n", DESIGN, R, B_BOOT,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("departures: rel = %.2f, both inject sd = %.2f logit units; off-index gamma = %.3f; sd(eta) = %.2f

", REL, REL * sd_eta, GAM_OFF, sd_eta))
fmt <- summ
fmt$ece <- sprintf("%.1f%%", 100 * fmt$ece); fmt$mae_true <- sprintf("%.1f%%", 100 * fmt$mae_true)
fmt$slope <- sprintf("%.2f", fmt$slope); fmt$brier <- sprintf("%.4f", fmt$brier); fmt$auc <- sprintf("%.3f", fmt$auc)
fmt$hl_heldout_rej <- sprintf("%.0f%%", 100 * fmt$hl_heldout_rej)
fmt$sc_edge_rej <- ifelse(is.na(fmt$sc_edge_rej), "", sprintf("%.0f%%", 100 * fmt$sc_edge_rej))
fmt$sc_hl_rej   <- ifelse(is.na(fmt$sc_hl_rej),   "", sprintf("%.0f%%", 100 * fmt$sc_hl_rej))
print(fmt[, c("case", "recal", "ece", "slope", "brier", "mae_true", "auc",
              "hl_heldout_rej", "sc_edge_rej", "sc_hl_rej")], row.names = FALSE)
cat("\nColumns: ece = expected calibration error (10 bins) on a 20000 test sample;",
    "\n  slope = calibration slope on test; mae_true = mean |p_hat - TRUE p| (only a simulation can see this);",
    "\n  auc on test; hl_heldout_rej = HL on a fresh n=1000 sample rejects at 5% (a test of H0^cal);",
    "\n  sc_*_rej = shrink.gof on the TRAINING data rejects at 5% (the test of H0^str).\n")
