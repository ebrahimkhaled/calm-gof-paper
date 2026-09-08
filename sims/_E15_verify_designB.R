# Verification of the review finding that "design B" names TWO different coefficient scalings.
#   E4c3_final_size_glaucoma.R (Tables 2, 7, 8, 9, Fig 4/5, supplement): b0 * 2.31 / sqrt(sum(b0^2))
#       -> the EUCLIDEAN NORM of beta is 2.31
#   E10_heldout_validation.R, E12_referee_items.R, E13, E14a (Tables 3, 4, 5, 6 and the Section 4.2
#       scale-recovery figures): b0 * 2.31 / sqrt(b0' Sigma b0)
#       -> the STANDARD DEVIATION OF THE LINEAR PREDICTOR is 2.31
# If both hold, the two are different designs and the paper describes only the first.
for (p in c(40, 100, 160)) {
  Sig <- 0.7^abs(outer(1:p, 1:p, "-"))
  b   <- rep_len(c(.35, -.3, .25, -.2, .15), p)
  bB  <- b * 2.31 / sqrt(sum(b^2))                        # design B  (E4c3)
  bP  <- b * 2.31 / sqrt(drop(t(b) %*% Sig %*% b))        # design B' (E10/E12/E13/E14a)
  cat(sprintf("p = %3d | design B : ||beta|| = %.2f, sd(eta) = %.3f | design B': ||beta|| = %.2f, sd(eta) = %.2f | signal ratio %.2f\n",
              p, sqrt(sum(bB^2)), sqrt(drop(t(bB) %*% Sig %*% bB)),
              sqrt(sum(bP^2)), sqrt(drop(t(bP) %*% Sig %*% bP)),
              sqrt(drop(t(bP) %*% Sig %*% bP)) / sqrt(drop(t(bB) %*% Sig %*% bB))))
}

cat("\n--- the claimed size collision at kappa = 0.25, lambda = 416 ---\n")
a <- read.csv("../results/E4c3_size_summary.csv")     # design B
b <- read.csv("../results/E10_heldout_summary.csv")   # design B'
print(a[a$cell == "B_p100_lam416", ], row.names = FALSE)
print(b[b$cell %in% c("H6_k25_weak", "B_k25", "H9_k25_G20"), ], row.names = FALSE)

cat("\n--- #11: which trace normalises each fitted constant? ---\n")
d1 <- read.csv("../results/delta_law_cells.csv")
cat("delta_law_cells.csv columns:", paste(names(d1), collapse = ", "), "\n")
m <- lm(meanv2_tr ~ 0 + kappa, data = d1)
cat(sprintf("refit of meanv2_tr (normalised by tr(Omega_MLE)) on kappa: coef %.4f, uncentred R2 %.3f\n",
            coef(m), 1 - sum(resid(m)^2) / sum(d1$meanv2_tr^2)))
if (file.exists("../results/delta_law_edge_cells.csv")) {
  d2 <- read.csv("../results/delta_law_edge_cells.csv")
  cat("delta_law_edge_cells.csv columns:", paste(names(d2), collapse = ", "), "\n")
  if ("g_dec" %in% names(d2)) {
    m2 <- lm(g_dec ~ 0 + kappa, data = d2)
    cat(sprintf("refit of g_dec (normalised by tr(Sigma_dn)) on kappa:    coef %.4f, uncentred R2 %.3f\n",
                coef(m2), 1 - sum(resid(m2)^2) / sum(d2$g_dec^2)))
  }
  if ("g_edge3" %in% names(d2)) {
    m3 <- lm(g_edge3 ~ 0 + kappa, data = d2)
    cat(sprintf("refit of g_edge3 (normalised by tr(P_Z Sigma_dn)) on kappa: coef %.4f, uncentred R2 %.3f\n",
                coef(m3), 1 - sum(resid(m3)^2) / sum(d2$g_edge3^2)))
  }
}
