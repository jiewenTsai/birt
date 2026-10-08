# Example 4. Quantile effects from distributional models (RTMB, ELGM): SHASH residuals,
# moderated scales, quantiles(), quantile_test() and simultaneous bands
library(birt)
d <- sim_rtirt(N = 500, K = 10, beta = c(0.15, 0))
y <- paste0("y", 1:10); t <- paste0("t", 1:10)

# 1. speed: SHASH residual (skewness, tail weight) in the latent regression
reg <- rtirt_syntax(y, t, structure = "path", cov = c("x1", "x2"))     # speed ~ ability + x1 + x2
fn <- birt(reg, d)
fs <- birt(reg, d, family = list(speed = shash()), control = rtmb_control(k_hyper = 13))
# (with the default k_hyper = 9 a warning may say that much posterior mass is on the outermost
#  outer nodes; more nodes per hyperparameter remove it)
convergence(fs)                         # mass on the outermost nodes vs its Gaussian value
compare(normal = fn, shash = fs)        # log marginal likelihoods
e <- estimates(fs); e[e$lhs == "speed" & e$op %in% c("sd", "skew", "tail"), ]
q <- quantiles(fs, p = c(0.1, 0.5, 0.9)); q[q$target == "speed", ]   # equal at every p: shape only

# 2. log RT: cross-loadings, and ability moderating the log scale of each RT
m5 <- paste(rtirt_syntax(y, t, cross = "free"),
            sprintf("V(%s) ~ ability", paste(t, collapse = " + ")), sep = "\n")
f5 <- birt(m5, d)
q5 <- quantiles(f5, p = c(0.1, 0.5, 0.9))
q5[q5$predictor == "ability", ]          # est, posterior SD, pointwise interval, simultaneous band
quantile_test(f5, p = c(0.1, 0.5, 0.9))  # does the effect change with p? (Wald test of V(t) ~ ability)
plot(f5, predictors = "ability")         # effect of ability on each quantile of log RT

# 3. SHASH with skewness moderated by ability (about a minute)
f5s <- birt(m5, d, family = list(.continuous = shash(skew = "ability")))
compare(normal = f5, shash = f5s)
quantile_test(f5s, p = c(0.1, 0.5, 0.9))

# 4. JAGS counterpart: asymmetric Laplace (ALD) fits, one per quantile
if (requireNamespace("rjags", quietly = TRUE)) {
  qa <- birt(reg, d, family = list(speed = ald(c(0.1, 0.5, 0.9))), engine = "jags",
             n_chains = 2, n_iter = 4000, n_burn = 2000, fit_indices = FALSE, progress = FALSE)
  print(qa)
  fv <- birt(paste(reg, "V(speed) ~ ability + x1 + x2", sep = "\n"), d, family = list(speed = shash()),
             control = rtmb_control(method = "aghq"))
  plot(fv, targets = "speed", ald = qa)  # RTMB quantile curves with the ALD estimates
}
