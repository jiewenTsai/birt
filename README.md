# birt

**BIvariate Response Time Modeling with 'lavaan'-Style Syntax**

[中文說明（繁體）](README.zh-TW.md)

`birt` fits item response theory models with one or two latent variables, for example
*ability* and *speed* in a joint model of response accuracy and response times, from
lavaan-style syntax. Binary items, ordinal (Likert) items and continuous indicators such as
log response times can be mixed in one model.

By default `birt(model, data)` integrates the persons out by adaptive Gauss-Hermite quadrature
and approximates the posterior of the remaining parameters as an extended latent Gaussian model
(ELGM; Stringer, Brown & Stafford, 2023) in 'RTMB'. It returns posterior means, SDs, 95%
intervals and the log marginal likelihood. The same marginal likelihood can be maximized
instead (standard errors, AIC / BIC, score-based invariance tests), and an optional JAGS engine
samples the full posterior for spike-and-slab priors and quantile (asymmetric Laplace) models.

## Installation

```r
install.packages("birt")                       # once on CRAN
# install.packages("remotes")
remotes::install_github("jiewenTsai/birt")     # development version
```

Required packages ('RTMB', 'lavaan', ...) are installed with it. Two optional extras:

```r
# engine = "jags": install JAGS >= 4.3 from https://mcmc-jags.sourceforge.io/, then
install.packages("R2jags")
# dif_tree():
install.packages("partykit")
```

## Quick start

```r
library(birt)
d <- sim_rtirt(N = 500, K = 10)        # y1..y10 (0/1), t1..t10 (log RT), covariates x1, x2
items <- paste0("y", 1:10); times <- paste0("t", 1:10)

model <- rtirt_syntax(items, times)    # van der Linden (2007); or write the syntax by hand
cat(model)
#> ability =~ y1 + y2 + y3 + y4 + y5 + y6 +
#>          y7 + y8 + y9 + y10
#> speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5 + -1*t6 +
#>          -1*t7 + -1*t8 + -1*t9 + -1*t10
#> ability ~~ speed

fit <- birt(model, d)                  # ELGM (default), about 20 s
fit                                    # overview (abridged)
#> birt 0.1.0, RTMB engine: marginal likelihood (persons integrated out)
#>   500 persons; latent variables: ability, speed
#>   Approximate Bayes (ELGM; W = person latents): 42 parameters
#>     inner: AGHQ (k = 9) per person; outer: AGHQ, k = 9 on lsd2, atr, Gaussian in 40 directions (81 nodes)
#>   log marginal likelihood -6883.73; -2logL at the posterior mode 13448.7
#>   Reliability: ability 0.760, speed 0.861

summary(fit)          # item parameters (IRT a, b), lavaan-style parameter table with the prior of each parameter
scores(fit)           # EAP and posterior SD per person; id = "..." carries a person identifier
reliability(fit)      # EAP reliability of each latent variable
```

Real data usually come as raw seconds with a person id:

```r
fit <- birt(rtirt_syntax(items, times), raw, log_rt = TRUE, id = "IDSTUD")
```

## What you can do

### One syntax, three ways to fit it

```r
fit_b  <- birt(model, d)                                           # ELGM: posterior means, SDs, 95% intervals, log marginal likelihood
fit_ml <- birt(model, d, control = rtmb_control(method = "aghq"))  # maximum likelihood: SEs, AIC, BIC, casewise scores
fit_j  <- birt(model, d, engine = "jags")                          # MCMC: 4 chains, DIC / WAIC / LOOIC, ppc()
```

All three share `summary()`, `estimates()`, `coef()`, `scores()`, `reliability()`,
`fit_indices()`, `convergence()` and `equations()`. The printed header says whether a fit uses
the marginal likelihood (persons integrated out) or the conditional likelihood (persons sampled).

### Compare models

```r
fc <- birt(rtirt_syntax(items, times, cross = "free"), d)       # ability also loads on each log RT
fi <- birt(rtirt_syntax(items, times, structure = "none"), d)   # ability and speed independent
compare(correlated = fit, cross = fc, independent = fi)         # log marginal likelihoods, Bayes factors
compare(correlated = fit_ml, cross = birt(rtirt_syntax(items, times, cross = "free"), d,
                                          control = rtmb_control(method = "aghq")))   # -2logL, AIC, BIC
```

With ELGM fits `compare()` lists the log marginal likelihood of each model and the Bayes factor
of the best model against each. With maximum likelihood fits it lists AIC and BIC. In the
simulated data the even items have a true cross-loading of -0.3, so the cross-loading model wins
by a Bayes factor of about 10^65 over the correlated model, and by 354 AIC points.

### Quantile effects from one distributional model

Let ability moderate the log scale of each response time (`V(t) ~ ability`), or use
sinh-arcsinh residuals (`shash()`) whose skewness depends on ability. One fit then gives the
effect of ability on *every* conditional quantile of each log response time.

```r
m5 <- paste(rtirt_syntax(items, times, cross = "free"),
            sprintf("V(%s) ~ ability", paste(times, collapse = " + ")), sep = "\n")
f5 <- birt(m5, d)
quantiles(f5, p = c(0.1, 0.5, 0.9))         # est, SD, pointwise interval and simultaneous band per quantile
quantile_test(f5, p = c(0.1, 0.5, 0.9))     # does the effect change with p? (Wald test of the scale moderation)
plot(f5, predictors = "ability")            # one panel per item: effect against p with a 95% band

f5s <- birt(m5, d, family = list(.continuous = shash(skew = "ability")))   # skewness moderated by ability
compare(normal = f5, shash = f5s)
cond_reliability(f5)                        # reliability of ability given its value: accuracy, RT, both
rt_information(f5)                          # where the RT information lies; what a time limit keeps
```

For item `t2` (true cross-loading -0.3) the effect of ability is about -0.26 at p = 0.1 and
-0.31 at p = 0.9: able respondents are faster, slightly more so in the slow tail. The
simultaneous band shows at which quantiles the curve excludes 0, and `quantile_test()` reports
that the change over p is not significant in these homogeneous data.

### Measurement invariance from one fit

Score-based tests (Merkle & Zeileis, 2013) use the casewise scores of a maximum likelihood
fit, so DIF and item drift are tested without refitting per group.

```r
ml <- birt(model, d, control = rtmb_control(method = "aghq"))
score_test(ml, d$x1 > 0, by_item = TRUE)    # DIF between two groups, item by item, Holm-adjusted
score_test(ml, rowSums(d[times]))           # drift of the item parameters along the total log time

z  <- data.frame(group = factor(d$x1 > 0), x2_high = factor(d$x2 > 0))
dif_tree(paste("ability =~", paste(items, collapse = " + ")), d, z)   # DIF tree (partykit)
```

Each test returns the statistic and p-value per parameter block (and per item with
`by_item = TRUE`); the tree splits the sample where the item parameters change most.

### Likert items moderated by a response time (MNLFA)

With `ordered = TRUE` the items follow a graded response model, and moderation statements
regress conditional moments on a person-level variable, as in moderated nonlinear factor
analysis: `E(item) ~ z` shifts thresholds (uniform DIF), `E(item) ~ z:f` moderates the
discrimination (nonuniform DIF), `E(f) ~ z` and `V(f) ~ z` the latent mean and variance.

```r
dh <- sim_hgrm(N = 600, J = 8)              # y1..y8 on a 1-4 scale; logT = log total response time
dh$logT <- as.numeric(scale(dh$logT))       # moderators are used as given, so standardize them
m <- '
  SA =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8
  E(y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8) ~ pa*logT:SA   # one common effect of time on every discrimination
  E(y1 + y2) ~ logT                                        # threshold shift in y1, y2; the other items are anchors
  E(SA) ~ logT                                             # latent mean (impact)
  V(SA) ~ logT                                             # latent variance
'
fo <- birt(m, dh, ordered = TRUE)
summary(fo)$moderation                      # each moderation effect with its 95% interval
summary(fo)$precision                       # precision of the scores by tercile of logT
ordinal_information(fo); plot(fo, type = "dif")
```

The output shows that fast respondents are measured less precisely: in the simulated data the
fastest third has a mean posterior SD of 0.58 and reliability .68, the slowest third 0.38 and
.86. `hgrm()` writes the same syntax from a few options:

```r
h <- hgrm("SA =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8", dh, moderators = "logT",
          a = "common", b = "free", anchor = paste0("y", 3:8), impact = c("mean", "var"))
cat(h$model)
```

### Spike-and-slab selection and quantile regression (JAGS)

```r
fs <- birt(rtirt_syntax(items, times, cross = "ssp"), d, engine = "jags")   # spike-and-slab on every cross-loading
summary(fs)$selection                                                      # p_incl, BF10, selected (p_incl > .5)

reg <- rtirt_syntax(items, times, structure = "path", cov = c("x1", "x2"))  # speed ~ ability + x1 + x2
fq  <- birt(reg, d, family = list(speed = ald(c(0.1, 0.5, 0.9))), engine = "jags")   # quantile regression of speed
fq                                                                         # one coefficient table per quantile
ppc(fs)                                                                    # posterior predictive checks
```

The selection table gives, for each cross-loading, the posterior probability that it is not
zero (`p_incl`) and the Bayes factor for inclusion (`BF10`). In the simulated data the odd items
(true value 0) get `p_incl` around .1 to .2 and are set to zero; the even items get `p_incl` = 1
with estimates near -0.3.

### Write up the model

```r
equations(fit)             # the model as equations, with priors; format = "latex" for an appendix
algorithm(fit)             # how the fit was computed, step by step
rtmb_code(ml)              # a stand-alone RTMB script that reproduces -2logL
jags_code(fs)              # the JAGS model file
```

## Syntax at a glance

| Syntax | Meaning |
|---|---|
| `f =~ y1 + y2` | Measurement loadings. A second `=~` line of the same factor adds cross-loadings |
| `f =~ -1*t1` | Fixed loading |
| `f =~ a*t1 + a*t2` | Equality constraint by label |
| `f =~ prior("ssp")*t1` | Spike-and-slab prior (JAGS) |
| `f =~ prior("dnorm(0, 100)")*t1` | Custom prior (JAGS precision parameterization) |
| `f1 ~~ f2` | Latent correlation (0 by default) |
| `f2 ~ f1 + x1` | Latent regression on a latent variable or covariates |
| `E(y) ~ z` | Threshold or intercept shift (uniform DIF) |
| `E(y) ~ z:f` | Loading moderated on the log scale (nonuniform DIF) |
| `V(t) ~ z` | Log residual SD of a continuous indicator |
| `E(f) ~ z`, `V(f) ~ z` | Latent mean and log SD |

Arguments of `birt()`: `ordered = TRUE` for Likert items, `family = list(speed = shash())` or
`list(.continuous = ald(0.5))` for residual distributions, `dp = dpriors(...)` for default
priors, `log_rt = TRUE` for raw seconds, `id = "..."` for person identifiers.

## Learn more

- [README.zh-TW.md](README.zh-TW.md): the full description in Traditional Chinese
- [birtExamples](https://github.com/jiewenTsai/birtExamples): tutorials and worked examples
- `vignette("birt")`, the tutorial `inst/tutorial/birt_tutorial.qmd` and the scripts
  `ex01`-`ex09` in `system.file("examples", package = "birt")`
- `?birt`, `?rtirt_syntax`, `?quantiles`, `?score_test`, `?hgrm`

## References

Bauer, D. J. (2017). A more general model for testing measurement invariance and differential
item functioning. *Psychological Methods, 22*, 507-526.

Jones, M. C., & Pewsey, A. (2009). Sinh-arcsinh distributions. *Biometrika, 96*, 761-780.

Merkle, E. C., & Zeileis, A. (2013). Tests of measurement invariance without subgroups: A
generalization of classical methods. *Psychometrika, 78*, 59-82.

Stringer, A., Brown, P., & Stafford, J. (2023). Fast, scalable approximations to posterior
distributions in extended latent Gaussian models. *Journal of Computational and Graphical
Statistics, 32*, 84-98.

van der Linden, W. J. (2007). A hierarchical framework for modeling speed and accuracy on test
items. *Psychometrika, 72*, 287-308.
