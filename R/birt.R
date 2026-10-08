#' Fit a Bayesian IRT / response time model from lavaan-style syntax
#'
#' By default the model is fitted in 'RTMB' by ELGM (extended latent Gaussian model; Stringer,
#' Brown & Stafford, 2023): approximate Bayesian inference on the marginal likelihood, the persons
#' integrated out by adaptive Gauss-Hermite quadrature and the remaining parameters by an outer
#' adaptive quadrature, with the priors of [dpriors()]. `rtmb_control(method = "aghq")` gives
#' maximum likelihood on the same marginal likelihood (needed by the score-based tools:
#' [score_test()], [quantile_score()], [dif_tree()]); `engine = "jags"` samples the conditional
#' likelihood (persons sampled with the parameters) with [R2jags::jags.parallel()], for
#' spike-and-slab priors and [ald()] quantile models.
#'
#' @param model Model syntax (see Details).
#' @param data A data frame with the indicators and covariates.
#' @param binary Names of the binary (0/1) indicators. By default every indicator whose
#'   observed values are all 0 or 1; all other indicators are continuous.
#' @param family Named list of residual distributions, e.g. `list(speed = ald(0.25))`;
#'   see [ald()]. Default: normal for every latent variable and continuous indicator.
#' @param n_chains,n_iter,n_burn,thin MCMC settings of `engine = "jags"`: chains (run in parallel, on at most
#'   `getOption("birt.cores")` processes; default the number of physical cores),
#'   iterations per chain, burn-in, thinning (default: 1000 kept draws per chain).
#' @param seed Random seed.
#' @param file Path for the generated JAGS model file (default: a temporary file).
#' @param label Optional label shown in the output.
#' @param fit_indices Compute the marginal-likelihood criteria after sampling.
#' @param progress Show a progress bar (runs the sampler in a background R session with
#'   'callr').
#' @param engine `"rtmb"` (default): the marginal likelihood (the persons integrated out) in
#'   'RTMB', by ELGM (default, see [rtmb_control()]: per-person AGHQ inside, adaptive quadrature
#'   over the hyperparameters outside; the same priors as JAGS, `dp` and `prior("...")`, no
#'   spike-and-slab) or by maximum likelihood (`control = rtmb_control(method = "aghq")`, which
#'   ignores the priors). `"jags"`: the conditional likelihood, Bayesian (the persons sampled
#'   together with the parameters). Instead of [ald()] the RTMB engine offers [shash()] residuals and
#'   moderated scales (`V(y) ~ z` in the syntax), and [quantiles()] gives the quantile effects at
#'   every level from one fit. Returns a `birt_rtmb` object.
#' @param control RTMB settings, [rtmb_control()].
#' @param log_rt `TRUE`: the continuous indicators are response times in seconds and are
#'   log-transformed before fitting (under the same names). A vector of column names
#'   log-transforms only those (e.g. when some continuous indicators are not times).
#' @param dp Default priors, [dpriors()] (engine `"jags"`, and engine `"rtmb"` with
#'   `method = "elgm"`); the prior of every parameter is shown in `summary()`.
#' @param ssp_p Prior inclusion probability of the spike-and-slab parameters: `NULL`
#'   (default) learns it (`p_incl ~ beta(1, 1)`, shared within loadings, within
#'   regressions and within the moderation effects of one kind and moderator); `BF10` then
#'   uses the posterior mean of the learned inclusion probability as the prior odds
#'   (`prior_incl` in `summary(fit)$selection`). A number fixes it, e.g. `0.5` for Bayes
#'   factors with prior odds 1, or `0.2` when few nonzero effects are expected.
#' @param id Person identifiers for [scores()]: a column name of `data` or a vector
#'   (default `1..N`).
#' @param ordered Ordinal (Likert) indicators: `TRUE` (every indicator with integer values)
#'   or their names. They can be combined with binary and continuous indicators (e.g.
#'   response times) and two latent variables; see "Ordinal indicators" and "Moderation".
#' @param itemtype Model of the ordinal indicators: `"grm"` (graded response, cumulative logits),
#'   `"gpcm"` (generalized partial credit, adjacent-category logits) or `"tppcm"` (two-parameter
#'   partial credit: one discrimination per step); one value, or a named vector per item. `"gpcm"`
#'   and `"tppcm"` are experimental: with moderation statements they are not yet checked by
#'   simulation-based calibration.
#' @param censor Censored continuous indicators (e.g. time limits): a named list with the limits
#'   of each, `list(t1 = c(upper = 60), t2 = c(lower = 1, upper = 60))`. A value at or above
#'   `upper` is right censored (only "at least `upper`" is known; it contributes the log survival
#'   function, log P(y >= upper | f), to the likelihood) and a value at or below `lower` is left
#'   censored (log P(y <= lower | f)); other values contribute their density, as in the usual
#'   likelihood of censored data (Meeker & Escobar, 1998; for time limits in IRT, Lee & Ying,
#'   2015). The limits are on the scale of the columns of `data` as given: seconds for the columns
#'   log-transformed by `log_rt`, whose limits are log-transformed with them; on the data's own
#'   scale otherwise (e.g. `log(60)` when `data` already holds log times). Normal and [shash()]
#'   residuals; `engine = "rtmb"` (ELGM and maximum likelihood). Values recorded beyond a limit are
#'   treated as censored at the limit. See "Censoring" in Details.
#'
#' @details
#' **Syntax.** At most two latent variables. Operators:
#' \tabular{ll}{
#' `f =~ y1 + y2` \tab loadings (the first `=~` line of a factor: measurement loadings,
#'   positive; later lines: cross-loadings) \cr
#' `f =~ -1*t1` \tab fixed loading \cr
#' `f =~ r*t1 + r*t2` \tab equality by label \cr
#' `f =~ prior("ssp")*t1` \tab spike-and-slab prior \cr
#' `f =~ prior("dnorm(0, 100)")*t1` \tab any JAGS prior (precision parameterization) \cr
#' `f1 ~~ f2` \tab factor correlation (default 0) \cr
#' `f ~~ 1*f` \tab fixed (residual) variance \cr
#' `f ~ x1 + prior("ssp")*x2` \tab latent regression on covariates \cr
#' `f2 ~ f1` \tab path between the latent variables \cr
#' }
#' A factor with at least one free loading has (residual) variance 1; a factor whose
#' loadings are all fixed (e.g. speed with loadings -1) has a free variance. Binary
#' indicators use a logit link, continuous indicators a normal (or ALD) distribution;
#' every indicator has an intercept. Covariates enter as given (standardize them first)
#' and may not be missing; indicators may be missing (at random).
#'
#' **Conditional and marginal likelihood.** The two engines treat the person latents
#' \eqn{\theta_j} differently:
#' \tabular{lll}{
#'   \tab conditional likelihood \eqn{p(y \mid \theta, \psi)} \tab marginal likelihood \eqn{\int p(y \mid \theta, \psi) p(\theta \mid \psi) d\theta} \cr
#' Bayesian \tab `engine = "jags"`: \eqn{\theta} sampled with the parameters \eqn{\psi} \tab `rtmb_control(method = "elgm")` \cr
#' maximum likelihood \tab (not offered: joint maximum likelihood is inconsistent) \tab `engine = "rtmb"` (AGHQ or Laplace) \cr
#' }
#' Both use the same model syntax and [equations()]; [algorithm()] shows the computation.
#' Person scores are posterior means and SDs of the sampled \eqn{\theta_j} (JAGS) or EAPs and
#' posterior SDs given the estimates (RTMB; under ELGM averaged over the outer nodes). The
#' fit indices of JAGS fits (DIC, WAIC, PSIS-LOO) use the marginal likelihood of each draw.
#' ("Conditional" here is about the person latents; [cond_reliability()] is reliability given
#' \eqn{\theta}, another use of the word.)
#'
#' **Identification with a path.** With unit speed loadings, a common cross-loading of
#' `theta` on the response times equals the path `speed ~ theta`, so a linear path
#' together with free cross-loadings on every response time is refused.
#'
#' **Ordinal indicators** (`ordered`). Categories recoded to 1..K per item; the loading is the
#' discrimination (positive: an ordinal item loads on one latent variable, on its first line).
#' `itemtype = "grm"`: cumulative logits \eqn{P(y_{ij} > k) = \mathrm{logit}^{-1}\{a_{ij}(f_j -
#' b_{ik} - \delta_{ij})\}}, thresholds \eqn{b_{i1} < b_{i2} < \dots} (Samejima, 1969);
#' `"gpcm"`: adjacent-category logits \eqn{\log P(y = k + 1) / P(y = k) = a_{ij}(f_j - b_{ik} -
#' \delta_{ij})} with unordered step difficulties (Muraki, 1992); `"tppcm"`: the same with one
#' discrimination per step, \eqn{a_{ik} = a_i r_{ik}}, \eqn{r_{i1} = 1} (Yu, 1991).
#'
#' **Moderation** (moderated nonlinear factor analysis; Bauer, 2017). A moderation statement
#' regresses a conditional moment on moderators z (observed columns of `data`, used as given):
#' `E()` the expectation, `V()` the variance.
#' \tabular{ll}{
#' `E(y1) ~ z` \tab shift of y1: categorical items on the latent scale,
#'   \eqn{\delta_{ij} = \beta_i z_j} (uniform DIF); continuous items their intercept \cr
#' `E(y1 + y2) ~ prior("ssp")*z` \tab spike-and-slab shifts (DIF screening; engine `"jags"`) \cr
#' `E(y1) ~ z:f` \tab loading of y1 on f, log scale, \eqn{\lambda_{ij} = \lambda_i e^{\alpha_i z_j}} (nonuniform DIF) \cr
#' `E(y1 + y2 + y3) ~ c1*z:f` \tab one common effect (equal labels), e.g. person precision \cr
#' `V(t1) ~ z` \tab log residual SD of a continuous indicator (z may be a latent variable, engine `"rtmb"`) \cr
#' `E(f) ~ z` (or `f ~ z`) \tab latent mean (impact) \cr
#' `V(f) ~ z` \tab latent SD (log scale) \cr
#' }
#' Indicators without a term on z (or with `0*z`) are anchors. `E(f) ~ z` together with
#' shifts of every indicator of f on z (free or equal) is not identified and is refused;
#' spike-and-slab shifts or at least one anchor identify it. [hgrm()] writes the usual
#' moderated graded response model in this syntax and fits it with `engine = "rtmb"`.
#'
#' **Default priors.** Measurement loadings hierarchical lognormal per latent variable
#' (`log lambda ~ N(mu_f, 1)`, `mu_f ~ N(0, 10)`; see [dpriors()]); cross-loadings and
#' regression coefficients `dnorm(0, 1)`; intercepts `dnorm(0, 1/9)` (binary) and
#' `dnorm(0, 1/100)` (continuous); residual SDs, factor SDs and ALD scales
#' half-t(3); correlation uniform(-1, 1); change them with `dp = dpriors(...)`. Spike-and-slab: `incl ~ dbern(p_incl)`,
#' `p_incl ~ beta(1, 1)` (or fixed, `ssp_p`), slab `N(0, sigma_slab^2)`, `sigma_slab ~ N+(0, 1)`, shared
#' within loadings, within regression coefficients and within the moderation effects of one
#' kind (shifts `E(y) ~ z` or log loadings `E(y) ~ z:f`) and moderator.
#'
#' **Censoring** (`censor`). With a time limit c the recorded time of a person who ran out of
#' time is not their response time; treating it as observed biases the intercepts, the residual
#' SDs and the speed variance downwards. `censor` replaces the density of such a value by the
#' probability of being beyond the limit, \eqn{S(c \mid f) = 1 - \Phi(T^{-1}((c - \mu_{ij}) /
#' \sigma_{ij}))} (\eqn{T} the identity for normal residuals, the standardized sinh-arcsinh
#' transform for [shash()]), inside the integral over the person latents. The number of censored
#' values is shown in the header; [cond_reliability()] of a censored fit gives the information
#' under its limits, and [rt_information()] the information kept under other limits.
#'
#' @references
#' Lee, Y.-H., & Ying, Z. (2015). A mixture cure-rate model for responses and response times in
#' time-limit tests. *Psychometrika, 80*(3), 748-775. \doi{10.1007/s11336-014-9419-8}
#'
#' Meeker, W. Q., & Escobar, L. A. (1998). *Statistical methods for reliability data*. Wiley.
#'
#' @return An object of class `birt` (a list), or `birt_quantile` when a vector of
#'   quantiles is given in `family`. See [birt-class] for its contents and accessors.
#' @examples
#' d <- sim_rtirt(N = 300, K = 5)
#' model <- rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5), cross = "ssp")
#' cat(model)
#' \donttest{
#' # JAGS engine (needs the JAGS library)
#' if (requireNamespace("rjags", quietly = TRUE)) {
#'   fit <- birt(model, data = d, engine = "jags", n_chains = 2, n_iter = 2000, n_burn = 1000,
#'               progress = FALSE)
#'   summary(fit)
#' }
#' # RTMB engine: maximum likelihood, and approximate Bayes (ELGM)
#' if (requireNamespace("RTMB", quietly = TRUE)) {
#'   m0 <- rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5))
#'   # maximum likelihood
#'   f_ml <- birt(m0, d, control = rtmb_control(method = "aghq"), progress = FALSE)
#'   f_bayes <- birt(m0, d, progress = FALSE)          # ELGM (default)
#'   summary(f_bayes)
#' }
#' }
#' @export
birt <- function(model, data, binary = NULL, family = NULL, n_chains = 4, n_iter = 10000,
                 n_burn = 5000, thin = NULL, seed = 2026, file = NULL, label = NULL,
                 fit_indices = TRUE, progress = TRUE, engine = c("rtmb", "jags"), control = rtmb_control(),
                 log_rt = FALSE, id = NULL, dp = dpriors(), ssp_p = NULL, ordered = NULL, itemtype = "grm",
                 censor = NULL) {
  cl <- match.call()
  engine <- match.arg(engine)
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  if (is.character(id) && length(id) == 1 && nrow(data) > 1) {
    if (!id %in% names(data)) stop("id column '", id, "' not in `data`")
    id <- data[[id]]
  }
  if (!is.null(id) && length(id) != nrow(data)) stop("`id` must have one value per row of `data`")
  if (isTRUE(log_rt) && !is.null(ordered) && !isFALSE(ordered))
    stop("with ordinal indicators give log_rt as the column names of the times to log-transform", call. = FALSE)
  data <- log_columns(data, log_rt, model)
  # several quantiles: one fit per quantile
  if (!is.null(family) && is.list(family) && !inherits(family, "birt_family")) {
    lens <- vapply(family, function(f) if (inherits(f, "birt_family")) length(f$p) else 1L, 1L)
    if (sum(lens > 1) > 1) stop("a vector of quantiles is allowed for one entry of `family` only")
    if (any(lens > 1)) {
      if (engine == "rtmb") stop("engine = \"rtmb\" fits one model for all quantiles: use shash() and quantiles(fit, p)", call. = FALSE)
      i <- which(lens > 1); ps <- family[[i]]$p
      fits <- lapply(ps, function(p) {
        fam <- family; fam[[i]] <- ald(p)
        f <- if (!is.null(file)) sub("(\\.jags)?$", sprintf("_p%03d.jags", round(1000 * p)), file) else NULL
        birt(model, data, binary, fam, n_chains, n_iter, n_burn, thin, seed, f, engine = "jags",
             label = sprintf("%s%s p = %s", label %||% "", if (is.null(label)) "" else ",", format(p)),
             fit_indices = fit_indices, progress = progress, id = id, dp = dp, ssp_p = ssp_p)
      })
      return(structure(list(fits = fits, p = ps, target = names(family)[i], call = cl), class = "birt_quantile"))
    }
  }
  if (!inherits(dp, "birt_dpriors")) stop("`dp` must come from dpriors()")
  if (!is.null(ssp_p) && !(is.numeric(ssp_p) && length(ssp_p) == 1 && ssp_p > 0 && ssp_p < 1)) stop("ssp_p must be NULL or a number in (0, 1)")
  spec <- build_spec(model, data, binary, family, ordered, itemtype)
  spec$log_rt <- attr(data, "birt_log_rt"); spec$dp <- dp; spec$ssp_p <- ssp_p
  if (!is.null(censor) && engine == "jags")
    stop("censor = needs engine = \"rtmb\" (the censored likelihood is not written into the JAGS model); fit with the default engine, by ELGM or rtmb_control(method = \"aghq\")", call. = FALSE)
  spec$censor <- censor_setup(censor, spec, data, spec$log_rt)
  if (!is.null(spec$ident)) {
    if (engine == "rtmb" && control$method != "elgm") stop(spec$ident$msg, call. = FALSE)
    warning(spec$ident$msg, "; here they are identified only through the priors", call. = FALSE)
  }
  if (engine == "rtmb") {
    mc <- intersect(names(cl), c("n_chains", "n_iter", "n_burn", "thin", "seed", "file", "fit_indices", "ssp_p"))
    if (length(mc)) message("engine = \"rtmb\" does not use ", paste(mc, collapse = ", "), " (settings of engine = \"jags\"); see rtmb_control()")
    return(birt_rtmb(cl, model, data, spec, control, label, id, progress))
  }
  if (!requireNamespace("R2jags", quietly = TRUE))
    stop("engine = \"jags\" needs JAGS (https://mcmc-jags.sourceforge.io) and the R2jags package: install.packages(\"R2jags\")", call. = FALSE)
  rt_only <- names(Filter(function(f) needs_rtmb(f, spec$factors), c(spec$family$factor, spec$family$ind)))
  if (length(rt_only)) stop("shash(), skew moderators and latent scale moderators need engine = \"rtmb\" (", paste(utils::head(rt_only, 3), collapse = ", "), ")", call. = FALSE)
  cg <- build_code(spec)
  if (is.null(file)) file <- tempfile(fileext = ".jags")
  writeLines(cg$code, file)
  dv <- unique(c(spec$binary, spec$cont, spec$covs, spec$modvars))
  jd <- c(lapply(stats::setNames(dv, dv), function(v) as.numeric(data[[v]])), N = nrow(data))
  for (v in spec$ordinal) jd[[v]] <- as.numeric(spec$ord$Y[, v])                   # categories 1..K
  if (is.null(thin)) thin <- max(1, floor((n_iter - n_burn) / 1000))

  est <- 20 + 5.5e-7 * jd$N * length(c(spec$binary, spec$cont)) * n_iter *
    (1 + sum(vapply(spec$family$ind, is_q, logical(1))) / max(1, length(spec$cont)))
  f <- run_bg("run_jags", list(file = file, data = jd, monitor = cg$monitor,
                               inits = make_inits(spec, data, cg, seed), n_chains = n_chains,
                               n_iter = n_iter, n_burn = n_burn, thin = thin, seed = seed),
              est_secs = est, label = if (is.null(label)) "Sampling" else sprintf("Sampling [%s]", label),
              show = progress)
  fit <- structure(list(call = cl, model = model, label = label, spec = spec, code = cg$code, file = file, likelihood = "conditional",
                        partab = cg$partab, data = jd, draws = f$draws, secs = f$secs,
                        mcmc = list(n_chains = n_chains, n_iter = n_iter, n_burn = n_burn, thin = thin, seed = seed),
                        id = id),
                   class = c("birt", "birt_fit"))
  fit <- add_derived(fit)
  if (fit_indices) fit <- add_fit_indices(fit, progress = progress)
  fit
}

# log-transform response time columns (names, or TRUE: every non-binary indicator of the model)
log_columns <- function(data, log_rt, model) {
  if (is.null(log_rt) || isFALSE(log_rt)) return(data)
  if (!isTRUE(log_rt) && !is.character(log_rt)) stop("log_rt must be TRUE / FALSE or column names", call. = FALSE)
  if (isTRUE(log_rt)) {
    vars <- parse_model(model)$ovs
    log_rt <- vars[!vapply(vars, function(v) all(data[[v]] %in% c(0, 1, NA)), TRUE)]
  }
  miss <- setdiff(log_rt, names(data))
  if (length(miss)) stop("log_rt: columns not in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  bad <- log_rt[vapply(log_rt, function(v) !is.numeric(data[[v]]) || any(data[[v]] <= 0, na.rm = TRUE), TRUE)]
  if (length(bad)) stop("log_rt: response times must be positive seconds; check ", paste(bad, collapse = ", "), call. = FALSE)
  for (v in log_rt) data[[v]] <- log(data[[v]])
  attr(data, "birt_log_rt") <- c(attr(data, "birt_log_rt"), log_rt)
  data
}

# number of parallel JAGS processes: option birt.cores (default: physical cores), at most 2
# under R CMD check (_R_CHECK_LIMIT_CORES_)
birt_cores <- function() {
  n <- getOption("birt.cores", parallel::detectCores(logical = FALSE))
  if (is.na(n) || n < 1) n <- 1
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(lim) && !identical(tolower(lim), "false")) n <- min(n, 2)
  as.integer(n)
}

run_jags <- function(file, data, monitor, inits = NULL, n_chains = 4, n_iter = 10000,
                     n_burn = 5000, thin = 5, seed = 2026) {
  # jags.parallel calls a no-argument inits function in each worker; each draws its own
  # chain id (different jitter). R2jags seeds the chains.
  ini <- if (is.null(inits)) NULL else {
    f0 <- inits
    function() f0(sample.int(1e6, 1))
  }
  t0 <- Sys.time()
  # partly defined arrays (thresholds of items with fewer categories, moderation effects of
  # some items only) are fine: R2jags' summary warns about them, the draws are complete
  fit <- withCallingHandlers(
    R2jags::jags.parallel(data = data, inits = ini, parameters.to.save = monitor,
                          model.file = file, n.chains = n_chains, n.iter = n_iter,
                          n.burnin = n_burn, n.thin = thin, n.cluster = min(n_chains, birt_cores()),
                          DIC = FALSE, jags.seed = seed, envir = new.env()),   # R2jags writes the data to `envir`: not the user's workspace
    warning = function(w) if (grepl("error/missing in parameter", conditionMessage(w))) invokeRestart("muffleWarning"))
  draws <- coda::as.mcmc(fit)
  draws <- draws[, setdiff(coda::varnames(draws), "deviance"), drop = FALSE]
  list(draws = draws, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# posterior summaries of the model parameters, person scores, reliability, item tables
add_derived <- function(fit) {
  M <- as.matrix(fit$draws); sp <- fit$spec
  person <- grepl(sprintf("^(%s)\\[", paste(sp$factors, collapse = "|")), colnames(M))
  pn <- colnames(M)[!person]
  P <- M[, pn, drop = FALSE]
  free <- apply(P, 2, stats::sd) > 0 & !grepl("^incl_", pn)       # no Rhat for 0/1 indicators
  rh <- ess <- stats::setNames(rep(NA_real_, length(pn)), pn)
  if (any(free)) {
    sub <- fit$draws[, pn[free], drop = FALSE]
    if (coda::nchain(sub) > 1) rh[free] <- coda::gelman.diag(sub, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]   # one chain: no R-hat
    ess[free] <- round(coda::effectiveSize(sub))
  }
  fit$summary_nodes <- data.frame(node = pn, mean = colMeans(P), sd = apply(P, 2, stats::sd),
                                  q025 = apply(P, 2, stats::quantile, 0.025), q975 = apply(P, 2, stats::quantile, 0.975),
                                  rhat = rh, ess = as.integer(ess), row.names = NULL)
  fit$estimates <- make_estimates(fit)
  fit$estimates <- ald_adjust(fit)
  N <- fit$data$N
  sc <- data.frame(id = fit$id %||% seq_len(N))
  rel <- numeric()
  for (f in sp$factors) {
    D <- M[, sprintf("%s[%d]", f, seq_len(N)), drop = FALSE]
    sc[[f]] <- colMeans(D); sc[[paste0(f, "_psd")]] <- apply(D, 2, stats::sd)
    rel[f] <- emp_rel(D)
  }
  fit$scores <- sc; fit$reliability <- rel
  fit$items <- make_items(fit)
  mod_outputs(fit)
}

#' Compute and store the marginal-likelihood fit indices
#'
#' Integrates the latent variables out (see Details of [birt()]) to get the pointwise
#' log-likelihood, then DIC, WAIC and PSIS-LOO. Called by [birt()] unless
#' `fit_indices = FALSE`.
#' @param fit A `birt` object.
#' @param n_draws Number of posterior draws used (reduced automatically when ALD adds
#'   quadrature layers).
#' @param Q Gauss-Hermite nodes.
#' @param progress Show progress.
#' @return `fit` with `fit_indices`, `loglik`, `loo` and `waic` added.
#' @export
add_fit_indices <- function(fit, n_draws = 1000, Q = 41, progress = TRUE) {
  nd <- n_draws_ll(fit$spec, n_draws)
  P <- prep_pars(fit, nd)
  P$f1_mom <- f1_moments(P); P$M <- NULL                     # keep only the centres of the adaptive rule
  chain <- P$chain
  res <- tryCatch(run_bg("fit_indices_compute", list(P = P, obs = obs_matrices(fit), chain = chain, Q = Q),
                         label = "Fit indices", show = progress),
                  error = function(e) { message("fit indices not computed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) {
    fit$fit_indices <- res$table; fit$loglik <- res$loglik; fit$loo <- res$loo; fit$waic <- res$waic
    attr(fit$loglik, "chain") <- chain
  }
  fit
}

# Background run with a progress display: fun (a birt function, by name) runs in a
# separate R session (callr::r_bg) and a spinner, elapsed time and a rough bar are shown
# (an estimate; it stops at 99% until the job ends). Without callr, with show = FALSE, or
# when birt is not installed (e.g. devtools::load_all()), fun runs directly.
run_bg <- function(fname, args, est_secs = NA, label = "Sampling", show = TRUE, width = 30) {
  fun <- get(fname, envir = asNamespace("birt"))
  installed <- length(find.package("birt", lib.loc = .libPaths(), quiet = TRUE)) > 0 &&
    !(requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("birt"))
  if (!show || !requireNamespace("callr", quietly = TRUE) || !installed) {
    if (show) message(label, " ...")
    return(do.call(fun, args))
  }
  job <- callr::r_bg(function(fname, args) {
    fun <- utils::getFromNamespace(fname, "birt"); do.call(fun, args)
  }, args = list(fname = fname, args = args), supervise = TRUE)
  t0 <- Sys.time(); spin <- c("|", "/", "-", "\\"); i <- 0
  fmt <- function(s) sprintf("%d:%02d", as.integer(s) %/% 60, as.integer(s) %% 60)
  repeat {
    job$wait(500)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!job$is_alive()) break
    i <- i + 1
    if (is.na(est_secs)) cat(sprintf("\r%s %s  elapsed %s   ", label, spin[i %% 4 + 1], fmt(el)))
    else {
      p <- min(el / est_secs, 0.99); n <- round(p * width)
      cat(sprintf("\r%s %s [%s%s] %3d%%  elapsed %s  ~remaining %s   ", label, spin[i %% 4 + 1],
                  strrep("=", n), strrep(" ", width - n), round(100 * p), fmt(el),
                  if (el < est_secs) fmt(est_secs - el) else "soon"))
    }
    utils::flush.console()
  }
  cat(sprintf("\r%s done [%s] 100%%  elapsed %s%s\n", label, strrep("=", width), fmt(el), strrep(" ", 20)))
  job$get_result()
}
