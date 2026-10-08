# Example 3. Latent regression: covariates of speed and ability, and a path ability -> speed (ELGM)
library(birt)
d <- sim_rtirt(N = 500, K = 10, beta = c(0.15, 0))      # x1 affects speed, x2 does not
y <- paste0("y", 1:10); t <- paste0("t", 1:10)
model <- rtirt_syntax(y, t, structure = "path", cov = c("x1", "x2"), cov_on = "both")
cat(model)                     # speed ~ ability + x1 + x2; ability ~ x1 + x2 (no cross-loadings, so the path is identified)
fit <- birt(model, d)          # ELGM
e <- estimates(fit)
e[e$op == "~", c("lhs", "op", "rhs", "est", "sd", "q025", "q975", "prior")]

# does x2 matter? Bayes factor of the model without x2
fit0 <- birt(rtirt_syntax(y, t, structure = "path", cov = "x1", cov_on = "both"), d)
compare(with_x2 = fit, without_x2 = fit0)

# shrinkage instead of selection: a small-variance prior on every regression coefficient
fit_s <- birt(model, d, dp = dpriors(beta = "dnorm(0, 25)"))     # SD 0.2
e_s <- estimates(fit_s)
cbind(e[e$op == "~", c("lhs", "rhs")], default = e$est[e$op == "~"], sd0.2 = e_s$est[e_s$op == "~"])

# spike-and-slab selection (prior("ssp")) needs engine = "jags"
if (requireNamespace("rjags", quietly = TRUE)) {
  m_ssp <- paste(rtirt_syntax(y, t), 'speed ~ prior("ssp")*x1 + prior("ssp")*x2', sep = "\n")
  f_ssp <- birt(m_ssp, d, engine = "jags", n_chains = 2, n_iter = 4000, n_burn = 2000, progress = FALSE)
  print(summary(f_ssp)$selection)
}
