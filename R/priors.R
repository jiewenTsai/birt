#' Default priors of a JAGS fit
#'
#' As `dpriors()` in blavaan: the default prior of each parameter class, as JAGS
#' distributions (normal priors in the **precision** parameterization, `dnorm(mean, 1 / sd^2)`).
#' A `prior("...")` modifier in the model syntax overrides the default for that parameter
#' (loadings and regression coefficients); `prior("ssp")` gives it a spike-and-slab prior.
#' The prior of every parameter is listed in the `prior` column of `summary()` and
#' [estimates()].
#' @param loading Measurement loadings (the first `=~` line of a factor; positive). The default
#'   `"hier"` is a hierarchical lognormal prior per latent variable, \eqn{\log\lambda_i \sim
#'   N(\mu_f, 1)}, \eqn{\mu_f \sim N(0, 10)}: the spread of `dlnorm(0, 1)`, with a centre that
#'   adapts to the units of the indicators. A fixed prior such as `"dlnorm(0, 1)"` (the earlier default)
#'   pulls all loadings of a latent variable towards 1; with many indicators on another scale
#'   (log response times: loadings about 0.1-0.3) this shrinks the latent variance and biases
#'   its correlation (with 170 items and 723 persons: .66 against .52).
#' @param cross Cross-loadings (later `=~` lines).
#' @param beta Regression coefficients (`~`).
#' @param intercept_binary,intercept Intercepts of binary (logit) and continuous indicators.
#' @param resid_sd Residual SD of continuous indicators (for [ald()]: of the implied SD).
#' @param factor_sd SD of a latent variable with a free variance.
#' @param cor Correlation of the latent variables (`f1 ~~ f2`).
#' @param slab_sd SD of the spike-and-slab slab (shared within loadings / within regressions).
#' @param moderation Log-scale moderation effects: `E(y) ~ z:f` (log loading), `V() ~ z` (log
#'   SD) and the `skew =` moderators of [shash()].
#' @param threshold,threshold_step Ordinal indicators: the first threshold (graded response;
#'   every step difficulty of a partial credit item), and the positive steps between the next
#'   thresholds of a graded response item (lognormal; Bayesian engines).
#' @param step_ratio Ratio of a later step discrimination to the first one (`itemtype = "tppcm"`).
#' @param skew,tail Engine `"rtmb"` with `method = "elgm"` only: the SHASH skewness and the
#'   SHASH tail weight (`dlnorm(0, 4)`: log tail weight with SD 0.5, centred at the normal).
#' @details The same priors are used by the approximate Bayesian (ELGM) fit of the RTMB engine,
#'   `birt(..., engine = "rtmb", control = rtmb_control(method = "elgm"))`, which supports
#'   `dnorm`, `dlnorm`, `dt`, `dunif`, `dgamma`, `dexp`, `dbeta`, `ddexp` and `dlogis`, with
#'   truncation `T(a, b)`.
#' @return A named list of prior strings (class `birt_dpriors`).
#' @examples
#' dpriors()
#' dpriors(cross = "dnorm(0, 100)")       # BSEM-style small-variance cross-loadings (SD 0.1)
#' dpriors(loading = "dlnorm(0, 1)")       # the earlier fixed loading prior
#' @export
dpriors <- function(loading = "hier", cross = "dnorm(0, 1)", beta = "dnorm(0, 1)",
                    intercept_binary = "dnorm(0, 1/9)", intercept = "dnorm(0, 1/100)",
                    resid_sd = "dt(0, 1, 3) T(0,)", factor_sd = "dt(0, 1, 3) T(0,)", cor = "dunif(-1, 1)",
                    slab_sd = "dnorm(0, 1) T(0,)",
                    moderation = "dnorm(0, 1)", threshold = "dnorm(0, 1/4)", threshold_step = "dlnorm(-0.5, 1)",
                    step_ratio = "dlnorm(0, 1)", skew = "dnorm(0, 1)", tail = "dlnorm(0, 4)") {
  p <- list(loading = loading, cross = cross, beta = beta, intercept_binary = intercept_binary, intercept = intercept,
            resid_sd = resid_sd, factor_sd = factor_sd, cor = cor, slab_sd = slab_sd,
            moderation = moderation, threshold = threshold, threshold_step = threshold_step, step_ratio = step_ratio, skew = skew, tail = tail)
  bad <- names(p)[!vapply(names(p), function(n) { x <- p[[n]]; is.character(x) && length(x) == 1 &&
    (grepl("^d[a-z]+\\(.*\\)", trimws(x)) || (n == "loading" && identical(trimws(x), "hier"))) }, TRUE)]
  if (length(bad)) stop("dpriors(): give JAGS distributions such as \"dnorm(0, 1)\" (loading: also \"hier\") for: ", paste(bad, collapse = ", "), call. = FALSE)
  structure(lapply(p, trimws), class = "birt_dpriors")
}

#' @export
print.birt_dpriors <- function(x, ...) {
  print_table(data.frame(parameter = names(x), prior = unlist(x)))
  invisible(x)
}

# ---- JAGS prior strings as log densities (RTMB engine, method = "elgm") -------------------
# "dnorm(0, 1/9) T(0,)" -> list(f = function(x) log density (works on RTMB advectors),
# lo, hi = support after truncation). JAGS parameterizations: normal / lognormal / t /
# double exponential / logistic with precision tau.
prior_density <- function(txt) {
  txt <- trimws(txt)
  trc <- regmatches(txt, regexec("\\s*T\\(([^()]*)\\)$", txt))[[1]]
  m <- regmatches(txt, regexec("^(d[a-z]+)\\(([^()]*)\\)", sub("\\s*T\\([^()]*\\)$", "", txt)))[[1]]
  if (!length(m)) stop("cannot read the prior \"", txt, "\"", call. = FALSE)
  m[4:5] <- if (length(trc)) trc[1:2] else c("", "")
  num <- function(a) {
    a <- trimws(a); if (a == "") return(NA_real_)
    v <- tryCatch(eval(str2lang(a), baseenv()), error = function(e) NULL)
    if (!is.numeric(v) || length(v) != 1) stop("prior \"", txt, "\": arguments must be numbers", call. = FALSE)
    v
  }
  a <- unname(vapply(strsplit(m[3], ",")[[1]], num, 0)); dn <- m[2]
  need <- c(dnorm = 2, dlnorm = 2, dt = 3, dunif = 2, dgamma = 2, dexp = 1, dbeta = 2, ddexp = 2, dlogis = 2)
  if (!dn %in% names(need)) stop("prior \"", txt, "\": ", dn, " is not supported (", paste(names(need), collapse = ", "), ")", call. = FALSE)
  if (length(a) != need[[dn]]) stop("prior \"", txt, "\": ", dn, " takes ", need[[dn]], " arguments", call. = FALSE)
  # support, log density and cdf (for the truncation constant)
  sup <- switch(dn, dlnorm = , dgamma = , dexp = c(0, Inf), dunif = a[1:2], dbeta = c(0, 1), c(-Inf, Inf))
  f <- switch(dn,
    dnorm = function(x) 0.5 * log(a[2] / (2 * pi)) - 0.5 * a[2] * (x - a[1])^2,
    dlnorm = function(x) 0.5 * log(a[2] / (2 * pi)) - 0.5 * a[2] * (log(x) - a[1])^2 - log(x),
    dt = function(x) lgamma((a[3] + 1) / 2) - lgamma(a[3] / 2) + 0.5 * log(a[2] / (a[3] * pi)) -
      (a[3] + 1) / 2 * log(1 + a[2] * (x - a[1])^2 / a[3]),
    dunif = function(x) 0 * x - log(a[2] - a[1]),
    dgamma = function(x) a[1] * log(a[2]) + (a[1] - 1) * log(x) - a[2] * x - lgamma(a[1]),
    dexp = function(x) log(a[1]) - a[1] * x,
    dbeta = function(x) (a[1] - 1) * log(x) + (a[2] - 1) * log(1 - x) - lbeta(a[1], a[2]),
    ddexp = function(x) log(a[2] / 2) - a[2] * sqrt((x - a[1])^2 + 1e-8),
    dlogis = function(x) { u <- a[2] * (x - a[1]); log(a[2]) - u - 2 * log(1 + exp(-u)) })
  cdf <- switch(dn,
    dnorm = function(q) stats::pnorm(q, a[1], 1 / sqrt(a[2])),
    dlnorm = function(q) stats::plnorm(q, a[1], 1 / sqrt(a[2])),
    dt = function(q) stats::pt((q - a[1]) * sqrt(a[2]), a[3]),
    dunif = function(q) stats::punif(q, a[1], a[2]), dgamma = function(q) stats::pgamma(q, a[1], a[2]),
    dexp = function(q) stats::pexp(q, a[1]), dbeta = function(q) stats::pbeta(q, a[1], a[2]),
    ddexp = function(q) ifelse(q < a[1], 0.5 * exp(a[2] * (q - a[1])), 1 - 0.5 * exp(-a[2] * (q - a[1]))),
    dlogis = function(q) stats::plogis(q, a[1], 1 / a[2]))
  lo <- sup[1]; hi <- sup[2]; cst <- 0
  if (nzchar(m[4])) {
    tr <- strsplit(paste0(m[5], " "), ",")[[1]]
    if (length(tr) != 2) stop("prior \"", txt, "\": truncation must be T(lower, upper)", call. = FALSE)
    tl <- num(tr[1]); th <- num(tr[2])
    if (!is.na(tl)) lo <- max(lo, tl)
    if (!is.na(th)) hi <- min(hi, th)
    cst <- -log(cdf(hi) - cdf(lo))
  }
  list(f = function(x) f(x) + cst, lo = lo, hi = hi, txt = txt)
}

# Hierarchical lognormal prior of the positive loadings of one latent variable ("hier"):
#   log lambda_i ~ N(mu_f, 1),  mu_f ~ N(0, HIER_V)
# The spread is that of the fixed dlnorm(0, 1); only the centre is estimated. A free tau would put
# an infinite density at equal loadings (tau -> 0) once tau is integrated out: harmless for MCMC,
# but the posterior mode of the ELGM fit falls into it. With mu integrated out the prior of
# u = log lambda is N(0, I + HIER_V 11'), in closed form.
HIER_V <- 10
hier_jags <- function(f) sprintf("  mu_load_%s ~ dnorm(0, %s)", f, num(1 / HIER_V))
hier_lp <- function(u) {
  n <- length(u); ub <- sum(u) / n
  -n / 2 * log(2 * pi) - 0.5 * log(1 + n * HIER_V) - 0.5 * (sum((u - ub)^2) + n * ub^2 / (1 + n * HIER_V))
}
