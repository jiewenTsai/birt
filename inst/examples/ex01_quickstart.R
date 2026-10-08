# Example 1. Quick start: joint model of accuracy and response times by ELGM (the default)
#
# birt(model, data) fits by ELGM (extended latent Gaussian model; Stringer, Brown & Stafford,
# 2023): approximate Bayes on the marginal likelihood. The persons are integrated out by
# per-person adaptive Gauss-Hermite quadrature (AGHQ), and the posterior of the other parameters
# by an outer adaptive quadrature, with the priors of dpriors().
library(birt)
d <- sim_rtirt(N = 500, K = 10)        # y1..y10 (0/1), t1..t10 (log RT), x1, x2
head(d)

y <- paste0("y", 1:10); t <- paste0("t", 1:10)
model <- rtirt_syntax(y, t)            # van der Linden (2007): ability ~~ speed
cat(model)

fit <- birt(model, d, label = "hierarchical")   # ELGM, about 20 s
fit                     # overview: log marginal likelihood, inner / outer quadrature
summary(fit)            # posterior means, SDs, 95% intervals and the prior of each parameter
head(scores(fit))       # EAP and posterior SD per person
reliability(fit)        # var(EAP) / (var(EAP) + mean PSD^2)
convergence(fit)        # optimizer, outer nodes, quadrature check

# real data: times in seconds and a person id
# fit <- birt(rtirt_syntax(items, times), raw, log_rt = TRUE, id = "IDSTUD")
