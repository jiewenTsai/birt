# birt 0.1.0

First release. birt fits item response and response time models written in lavaan-style
syntax, with persons j and items i (y_ij is the response of person j to item i).

## Engines

* **ELGM (default).** `birt(model, data)` fits in 'RTMB' (`engine = "rtmb"`) with
  `rtmb_control(method = "elgm")`: approximate Bayesian inference for an extended latent
  Gaussian model (Stringer, Brown & Stafford, 2023). The persons are integrated out by
  per-person adaptive Gauss-Hermite quadrature, the other parameters by an outer adaptive
  quadrature at the posterior mode, with the priors of `dpriors()`. It gives posterior means,
  SDs and intervals, EAP scores and the log marginal likelihood for `compare()` (Bayes
  factors). `convergence()` reports the outer grid, and a warning flags a posterior that is
  too wide for the grid along a hyperparameter. Checked: with a Gaussian outer layer the log
  marginal likelihood equals the Laplace value at the mode (to 1e-8). `hyper = "auto"` puts
  outer quadrature on at most 3 hyperparameters (SHASH shape first) and says which it treats as
  Gaussian: with a 4th direction (5^4 = 625 nodes) the outer step took 5 times longer for
  posterior means within 0.034 SD. Timing (N = 1000-2000, 20-40 items, two latent variables,
  normal or SHASH; 81 outer nodes): ELGM took 1.4-2.0 times the maximum likelihood fit.
* **Maximum likelihood.** `rtmb_control(method = "aghq")`: the same marginal likelihood
  maximized, with standard errors, AIC / BIC and the casewise scores used by `score_test()`,
  `quantile_score()` and `dif_tree()`. Checked against brute-force integration (to 1e-4) and
  against mirt's GPCM and nominal models (-2logL within 1e-5).
* **JAGS (optional).** `engine = "jags"` (R2jags and rjags in Suggests) samples the person
  latents with the parameters (conditional likelihood): spike-and-slab priors, `ald()`
  quantile models, DIC / WAIC / PSIS-LOO (marginal likelihood per draw by adaptive
  Gauss-Hermite quadrature; Liu & Pierce, 1994), posterior predictive checks (`ppc()`).
  Starting values come from the data, and each chain starts the scale of each latent variable
  at a different place so that R-hat can see slow mixing along the scale.
* Printed headers, `equations()` and `fit$likelihood` say whether a fit uses the conditional
  or the marginal likelihood (Merkle, Furr & Rabe-Hesketh, 2019).

## Syntax

* `birt()`: one or two latent variables, binary, ordinal and continuous (e.g. log response
  time) indicators, latent regression, factor correlation, labels for equality constraints,
  `prior()` modifiers and fixed values. `rtirt_syntax()` writes the usual response time IRT
  models (van der Linden, 2007; cross-loadings, path or correlation structure, covariates).
* Default priors `dpriors()`: a hierarchical lognormal prior for the measurement loadings of
  each latent variable (log lambda ~ N(mu_f, 1), mu_f ~ N(0, 10)), and priors for thresholds,
  steps, moderation effects, SDs and the correlation; the prior of every parameter is shown in
  `summary()`.
* One way to write each thing: scale moderation is `V(y) ~ z`, all continuous indicators are
  `family = list(.continuous = )`. Ambiguous or ignored specifications are errors (a label
  with `prior()` on one term, labels shared across parameter types, a label times a number).
* Real data: `log_rt = TRUE` for raw seconds, `id =` for person identifiers, 1/2 coding of
  binary items, and a warning when continuous indicators look like raw times.

## Ordinal indicators and MNLFA

* `ordered = TRUE` (or names) makes integer indicators with 3-10 categories ordinal, combined
  freely with binary and continuous indicators, by every engine. `itemtype =` chooses the
  graded response model (default), the generalized partial credit model (Muraki, 1992) or the
  two-parameter partial credit model (Yu, 1991).
* Moderated nonlinear factor analysis (Bauer, 2017) for every indicator type: `E(y) ~ z`
  (threshold or intercept shift), `E(y) ~ z:f` (log loading), `V(y) ~ z` (log residual SD),
  `E(f) ~ z` and `V(f) ~ z`, with priors, labels, anchors (`0*z`) and spike-and-slab
  (`prior("ssp")`, JAGS). Moderation of the partial credit models is experimental (not yet
  checked by simulation-based calibration).
* `hgrm()` writes the moderated graded response model from the options `a`, `b`, `impact` and
  `anchor` and fits it by ELGM (or maximum likelihood for `score_test()`).
* `fit$dif`, `fit$impact`, `fit$precision` (score precision by moderator tercile),
  `ordinal_information()` and `plot(fit, type = "dif" / "information")`.

## Distributional and quantile tools

* `shash()` residuals (sinh-arcsinh; Jones & Pewsey, 2009) for continuous indicators and
  latent variables, with skewness moderation `shash(skew = )`.
* `quantiles()`: quantile effects of the predictors at every level p, with a simultaneous 95%
  band over the levels (sup-t; Montiel Olea & Plagborg-Moller, 2019) and their joint
  covariance (posterior covariance with ELGM). Checked against numerical derivatives of the
  conditional quantiles (to 1e-6).
* `quantile_test()`: does an effect change with p? It does exactly when the predictor
  moderates the scale or skewness of the target, so the test is the Wald test of those
  parameters (Koenker & Bassett, 1982). Checked: 5 of 56 rejections at .05 under the null.
* `quantile_score()`: casewise influence functions of the quantile effects (Hampel, 1974),
  robust standard errors, case influence, conditional or unconditional on covariates
  (Graubard & Korn, 1999), and `score_test(quantile_score(fit, p), z)` for invariance of an
  effect. Checked: robust SEs agree with 100 bootstrap fits (0.056 / 0.049 vs 0.058 / 0.051).
* `quantile_information()`: which quantile regression of the response times says the most about
  the first latent variable: the information of the tau-th regression (the Godambe information of
  its check-loss estimating function, f(q)^2 beta(tau)^2 / (tau (1 - tau)); Godambe, 1960;
  Koenker, 2005), speed projected out, and its efficiency against the whole RT distribution.
  Checked: beta equals the cross-loading and the item values equal the closed form for normal
  residuals; for SHASH beta and f(q) match the quantile function (to 1e-5) and the information
  agrees with Monte Carlo of the dichotomized times (within 3%).
* `ald()` (JAGS): Bayesian quantile regression with the asymmetric Laplace working likelihood
  (Yu & Zhang, 2005; Kozumi & Kobayashi, 2011), measurement or latent level, with intervals
  adjusted as in Yang, Wang & He (2016).

## Reliability and information

* `cond_reliability()`: conditional reliability of the first latent variable,
  I / (I + 1 / var) (Nicewander, 2018), with speed as a nuisance, the information of the
  response-time location, scale and skewness, and delta-method standard errors. Checked
  against the 2PL information, a closed form and Monte Carlo. Ordinal items (GRM, GPCM,
  TPPCM) and moderated loadings and shifts are included at the moderator values `at` (a
  censored fit: the information under its limits, Escobar & Meeker, 1998). Checked
  against `ordinal_information()` (to 1e-8) and the numerical Fisher information of the
  likelihood code (to 1e-6; censored and SHASH items to 1e-5).
* `rt_information()` (non-core): where in the response-time distribution the information lies (bands of
  RT quantile levels; Zheng & Gastwirth, 2000) and how much is kept when times beyond a limit,
  or fast responses, are censored (Efron & Johnstone, 1990). Checked: the bands sum to the RT
  information (to 3e-7) and the censored information agrees with Monte Carlo (within 0.1%)
  and with the curvature of the censored likelihood of `censor =` (to 1e-5). Moderated
  loadings and shifts and ordinal items are included, as in `cond_reliability()`.
* `ordinal_information()`: test information of ordinal indicators at chosen moderator values.
* Empirical reliability of the EAP scores in every fit.

## Invariance tools

* `score_test()`: score-based tests of parameter invariance (Merkle & Zeileis, 2013) from one
  maximum likelihood fit, with casewise scores from `estfun()`: DIF along a numeric, ordinal
  or categorical moderator, per item with Holm adjustment, a message when the moderator is
  related to the scores (impact), and cluster-level moderators. Valid for the measurement
  parameters in a joint accuracy and response-time model.
* The score tools are maximum likelihood methods: on an ELGM fit `score_test()`,
  `quantile_score()` and `estfun()` (and `dif_tree()` with an ELGM `control`) stop with the
  instruction to refit with `rtmb_control(method = "aghq")`. At the posterior mode the scores
  sum to minus the gradient of the log prior and MAP versions need centred scores with
  covariance-dependent error rates (Debelak, Pawel, Strobl & Merkle, 2022); the Bayesian
  route for one effect is `compare()` of ELGM fits with and without it (Bayes factor).
* `dif_tree()`: DIF trees by model-based recursive partitioning (Zeileis, Hothorn & Hornik,
  2008; Strobl, Kopf & Zeileis, 2015) for binary and ordinal items, with latent mean
  differences modelled in each node (`impact = TRUE`). Numeric partitioning variables are cut
  into `nbins = 10` quantile groups and tested as ordinal (Merkle, Fan & Zeileis, 2014):
  a 2PL tree (N = 600) with a numeric variable took 26 s instead of 15 minutes for the
  exhaustive search (`nbins = NULL`), and split at 6.01 for a true cut point at 6 (exhaustive:
  6.05).
* Spike-and-slab screening (JAGS) of loadings, regressions and moderation effects, with
  `p_incl`, the learned prior inclusion probability and `BF10` in `summary(fit)$selection`.

## Checks

* `check_sbc()`: simulation-based calibration of the JAGS engine (Talts et al., 2018), with
  an R simulator written independently of the generated JAGS code.
* `convergence()`: R-hat and ESS (JAGS); optimizer, gradient, quadrature accuracy and the
  outer grid (RTMB).

## Output helpers

* `summary()` returns the tables of a fit as one object (`$fit`, `$items`, `$se`,
  `$parameters` with the prior of each parameter, `$moderation`, `$selection`, ...), printed
  in lavaan / blavaan sections; `as.data.frame()` gives the parameter table.
* `equations()` writes the model as equations with the priors (Unicode, or LaTeX with
  `format = "latex"`); `algorithm()` writes how the fit was computed, with the numbers of the
  fit (JAGS samplers, AGHQ rounds, ELGM outer quadrature).
* `jags_code()` and `rtmb_code()` export a script that reproduces the fit, from a fit or from
  syntax (with `ordered =`, `itemtype =`). `rtmb_code()` writes ordinal items (GRM, GPCM, TPPCM),
  `E()` moderation and censored indicators; checked: the
  script reproduces -2logL of the fit (to 1e-5) and its joint density equals birt's (to 1e-10).
  `scores()`,
  `estimates()`, `posterior_draws()`, `compare()`, `plot()`.
* Generics shared with the sister package birtRcpp.

## Non-core features

* Censored response times (time limits): `birt(..., censor = list(t1 = c(upper = 60)))`, with
  limits on the scale of the data columns (log-transformed with the times under `log_rt`). A
  value at or beyond a limit contributes log S(c | f) or log F(c | f) instead of the density
  (Meeker & Escobar, 1998; Lee & Ying, 2015); normal and SHASH residuals, ELGM and maximum
  likelihood (not JAGS). Checked: -logL equals brute-force integration (normal to 1e-13; SHASH:
  the likelihood code on a fine grid to 1e-8); in a time-limit simulation (N = 1000, 10 items,
  18% censored, 8 replications) the mean bias of the time intercepts, residual SDs and speed SD
  is -0.0003, 0.0003 and 0.004, against -0.066, -0.063 and -0.079 when the censored times are
  treated as observed.
