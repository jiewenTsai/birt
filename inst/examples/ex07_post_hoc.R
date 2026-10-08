# Example 7. Using a fitted model for further quantities (no refitting)
library(birt)
d <- sim_rtirt(N = 500, K = 10)
y <- paste0("y", 1:10); t <- paste0("t", 1:10)
model <- rtirt_syntax(y, t, cross = "free")     # item-specific cross-loadings of ability on log RT
fit <- birt(model, d)                           # ELGM

# 1. scores: EAP and posterior SD per person (mixed over the outer ELGM nodes)
sc <- scores(fit); head(sc)

# 2. reliability in subgroups from the EAPs and posterior SDs
rel <- function(m, psd) var(m) / (var(m) + mean(psd^2))
g <- d$x1 > 0
c(low = rel(sc$ability[!g], sc$ability_psd[!g]), high = rel(sc$ability[g], sc$ability_psd[g]))

# 3. conditional reliability of ability: accuracy, response times, both
cr <- cond_reliability(fit)
cr
plot(cr)

# 4. where the response times carry information about ability, and what a time limit keeps
ri <- rt_information(fit)
ri
plot(ri)

# 5. parameters for reporting, and the posterior covariance of the quantile effects
summary(fit)$items
q <- quantiles(fit, p = c(0.25, 0.75))
head(q[q$predictor == "ability", ])
dim(attr(q, "vcov"))

# 6. JAGS: draws of every person's latent variables (plausible values, any function of them)
if (requireNamespace("rjags", quietly = TRUE)) {
  jf <- birt(model, d, engine = "jags", n_chains = 2, n_iter = 4000, n_burn = 2000, progress = FALSE)
  M  <- posterior_draws(jf, persons = TRUE)
  D  <- M[, grep("^ability\\[", colnames(M))]        # draws x persons
  pv <- D[sample(nrow(D), 5), ]                       # 5 plausible values per person
  print(dim(pv))
  r <- apply(D, 1, function(th) cor(th, d$x2))        # correlation with an external variable
  print(quantile(r, c(.025, .5, .975)))
  print(head(order(-jf$loo$diagnostics$pareto_k)))    # most influential persons (PSIS-LOO)
}
