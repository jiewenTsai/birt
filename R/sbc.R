# Simulation-based calibration of the JAGS engine (Talts, Betancourt, Simpson, Vehtari & Gelman, 2018).
# The parameters are drawn from the priors and the data from the model by an R simulator written
# from the model's definition (not from the generated JAGS code), so a mistake in the translation of
# the syntax into JAGS shows up as non-uniform ranks of the true values among the posterior draws.

#' Simulation-based calibration of a JAGS model
#'
#' Checks that a model is fitted correctly by the JAGS engine. Each replication draws the
#' parameters from their priors and data from the model (persons, covariates, moderators and
#' the pattern of missing values as in `data`), fits the model with [birt()] and records the
#' rank of every true parameter value among the posterior draws. If the model, the priors and
#' the sampler agree, each rank is uniform (Talts et al., 2018). The simulator is written in R
#' from the definition of the model, independently of the JAGS code that [birt()] generates, so
#' the check also tests the translation of the syntax into JAGS.
#'
#' Ties (spike-and-slab parameters at 0) get randomized ranks. Replications whose data would
#' change the model (an ordinal category or a binary value that is never observed) are drawn
#' again. Asymmetric Laplace families are not simulated.
#' @param model,data,binary,family,ordered,itemtype,dp,ssp_p As in [birt()]; `data` gives the
#'   persons, covariates, moderators and missing values.
#' @param n_rep Number of replications (100 or more for a real check; each one is a JAGS fit).
#' @param n_draws Posterior draws used for each rank (equally spaced, which also thins).
#' @param n_chains,n_iter,n_burn MCMC settings of each fit.
#' @param bins Number of bins of the uniformity test.
#' @param seed Random seed.
#' @param progress Show progress.
#' @return A `birt_sbc` object: `ranks` (replications x parameters), `test` (per parameter: the
#'   chi-square uniformity statistic and p-value, and the mean rank relative to its expectation),
#'   `redraws` and `rhat` (the largest R-hat of each replication; ranks are only meaningful when the
#'   chains converged, Talts et al., 2018, Section 5). `plot()` shows the rank histograms.
#' @references Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018).
#'   Validating Bayesian inference algorithms with simulation-based calibration. arXiv:1804.06788.
#' @examples
#' \donttest{
#' if (requireNamespace("rjags", quietly = TRUE)) {
#'   d <- sim_rtirt(N = 100, K = 3)
#'   m <- rtirt_syntax(paste0("y", 1:3), paste0("t", 1:3))
#'   s <- check_sbc(m, d, n_rep = 20, n_iter = 1000, n_burn = 500, progress = FALSE)
#'   s
#' }
#' }
#' @export
check_sbc <- function(model, data, n_rep = 100, n_draws = 100, binary = NULL, family = NULL, ordered = NULL, itemtype = "grm",
                      dp = dpriors(), ssp_p = NULL, n_chains = 2, n_iter = 2000, n_burn = 1000, bins = 10, seed = 1,
                      progress = TRUE) {
  local_seed(seed)
  sp <- build_spec(model, data, binary, family, ordered, itemtype); sp$dp <- dp; sp$ssp_p <- ssp_p
  if (any(vapply(c(sp$family$factor, sp$family$ind), is_q, TRUE))) stop("check_sbc(): ald() families are not simulated", call. = FALSE)
  if (any(vapply(c(sp$family$factor, sp$family$ind), function(f) needs_rtmb(f, sp$factors), TRUE))) stop("check_sbc() is for models of the JAGS engine", call. = FALSE)
  pt <- build_code(sp)$partab
  nodes <- unique(stats::na.omit(pt$node))
  ranks <- matrix(NA_real_, n_rep, length(nodes), dimnames = list(NULL, nodes)); redraws <- 0; rhat <- rep(NA_real_, n_rep)
  for (r in seq_len(n_rep)) {
    repeat {
      sim <- sbc_simulate(sp, pt, data)
      ok <- tryCatch({ s2 <- build_spec(model, sim$data, binary, family, ordered, itemtype)
                       identical(s2$binary, sp$binary) && identical(s2$ordinal, sp$ordinal) && identical(s2$ord$K, sp$ord$K) }, error = function(e) FALSE)
      if (ok) break
      redraws <- redraws + 1
      if (redraws > 20 * n_rep) stop("check_sbc(): the simulated data rarely have every category; use more persons", call. = FALSE)
    }
    fit <- suppressWarnings(suppressMessages(birt(model, sim$data, engine = "jags", binary = binary, family = family, ordered = ordered, itemtype = itemtype,
                                                  dp = dp, ssp_p = ssp_p, n_chains = n_chains, n_iter = n_iter, n_burn = n_burn,
                                                  seed = seed + r, fit_indices = FALSE, progress = FALSE)))
    D <- as.matrix(fit$draws)
    if (n_chains > 1) rhat[r] <- max(coda::gelman.diag(fit$draws[, nodes], autoburnin = FALSE, multivariate = FALSE)$psrf[, 1], na.rm = TRUE)
    keep <- unique(round(seq(1, nrow(D), length.out = min(n_draws, nrow(D)))))
    for (n in nodes) {
      x <- D[keep, n]; t0 <- sim$truth[[n]]
      ranks[r, n] <- sum(x < t0) + floor(stats::runif(1) * (sum(x == t0) + 1))           # ties at random
    }
    if (progress) message(sprintf("check_sbc(): replication %d of %d", r, n_rep))
  }
  test <- do.call(rbind, lapply(nodes, function(n) {
    x <- ranks[, n]
    br <- seq(-0.5, n_draws + 0.5, length.out = bins + 1); o <- table(cut(x, br, include.lowest = TRUE)); e <- n_rep / bins
    chi <- sum((o - e)^2 / e)
    data.frame(parameter = sbc_label(pt, n), node = n, chisq = chi, p = stats::pchisq(chi, bins - 1, lower.tail = FALSE),
               mean_rank = mean(x) / n_draws, stringsAsFactors = FALSE)
  }))
  structure(list(ranks = ranks, test = test, n_rep = n_rep, n_draws = n_draws, bins = bins, redraws = redraws, rhat = rhat), class = "birt_sbc")
}

sbc_label <- function(pt, n) { i <- which(pt$node == n)[1]; if (pt$op[i] %in% c("|", "|a")) paste0(pt$lhs[i], pt$op[i], pt$rhs[i]) else paste(pt$lhs[i], pt$op[i], pt$rhs[i]) }

#' @rdname check_sbc
#' @param x A `birt_sbc` object.
#' @param ... Unused.
#' @export
print.birt_sbc <- function(x, ...) {
  t <- x$test; bad <- t[t$p < 0.01, ]
  cat(sprintf("Simulation-based calibration: %d replications, ranks among %d draws; %d parameters\n", x$n_rep, x$n_draws, nrow(t)))
  cat(sprintf("Uniformity (chi-square, %d bins): smallest p = %.3g; %d parameters with p < .01 (expected about %.1f by chance); Holm-adjusted smallest p = %.3g\n",
              x$bins, min(t$p), nrow(bad), 0.01 * nrow(t), min(1, min(t$p) * nrow(t))))
  if (x$redraws) cat(sprintf("(%d data sets drawn again because a category was not observed)\n", x$redraws))
  nb <- sum(x$rhat > 1.05, na.rm = TRUE)
  if (nb) cat(sprintf("Not converged (max R-hat > 1.05) in %d of %d replications (largest %.2f): non-uniform ranks then reflect the sampler, not the model code; use more iterations or a better identified model\n",
                      nb, x$n_rep, max(x$rhat, na.rm = TRUE)))
  print_table(t[order(t$p), c("parameter", "chisq", "p", "mean_rank")][seq_len(min(10, nrow(t))), ], 3)
  invisible(x)
}

#' @rdname check_sbc
#' @export
plot.birt_sbc <- function(x, ...) {
  t <- x$test[order(x$test$p), ][seq_len(min(12, nrow(x$test))), ]
  nc <- ceiling(sqrt(nrow(t))); op <- graphics::par(mfrow = c(ceiling(nrow(t) / nc), nc), mar = c(3, 3, 2, 1)); on.exit(graphics::par(op))
  for (i in seq_len(nrow(t))) {
    graphics::hist(x$ranks[, t$node[i]], breaks = seq(-0.5, x$n_draws + 0.5, length.out = x$bins + 1), main = sprintf("%s (p = %.2f)", t$parameter[i], t$p[i]),
                   xlab = "rank", ylab = "", cex.main = 0.8, ...)
    graphics::abline(h = x$n_rep / x$bins, lty = 2)
  }
  invisible(x)
}

# ---- draws from the JAGS prior strings (precision parameterization, truncation T(a, b)) --------
rprior <- function(txt, n = 1) {
  d <- prior_density(txt)
  m <- regmatches(txt, regexec("^(d[a-z]+)\\(([^()]*)\\)", sub("\\s*T\\([^()]*\\)$", "", trimws(txt))))[[1]]
  a <- vapply(strsplit(m[3], ",")[[1]], function(x) eval(str2lang(trimws(x)), baseenv()), 0)
  q <- switch(m[2],
    dnorm = function(u) stats::qnorm(u, a[1], 1 / sqrt(a[2])), dlnorm = function(u) stats::qlnorm(u, a[1], 1 / sqrt(a[2])),
    dt = function(u) a[1] + stats::qt(u, a[3]) / sqrt(a[2]), dunif = function(u) stats::qunif(u, a[1], a[2]),
    dgamma = function(u) stats::qgamma(u, a[1], a[2]), dexp = function(u) stats::qexp(u, a[1]),
    dbeta = function(u) stats::qbeta(u, a[1], a[2]), dlogis = function(u) stats::qlogis(u, a[1], 1 / a[2]),
    stop("check_sbc(): no sampler for the prior ", txt, call. = FALSE))
  cdf <- function(x) { if (is.infinite(x)) as.numeric(x > 0) else stats::uniroot(function(u) q(u) - x, c(1e-15, 1 - 1e-15), tol = 1e-12)$root }
  ul <- if (is.finite(d$lo) && d$lo > q(1e-15)) cdf(d$lo) else 0
  uh <- if (is.finite(d$hi) && d$hi < q(1 - 1e-15)) cdf(d$hi) else 1
  q(stats::runif(n, ul, uh))
}

# ---- one data set from the model: parameters from the priors, then persons and indicators -------
sbc_simulate <- function(sp, pt, data) {
  dp <- sp$dp; N <- nrow(data); truth <- list(); fs <- sp$factors
  setv <- function(node, v) { truth[[node]] <<- v; v }
  ld <- sp$ld
  # spike-and-slab hyperparameters per group
  grp <- unique(stats::na.omit(c(ifelse(!is.na(pt$incl) & pt$kind == "loading", "load", NA), ifelse(!is.na(pt$incl) & pt$kind == "beta", "reg", NA), pt$group)))
  pin <- list(); ssd <- list()
  for (g in grp) {
    pin[[g]] <- if (is.null(sp$ssp_p)) setv(sprintf("p_incl_%s", g), stats::rbeta(1, 1, 1)) else sp$ssp_p
    ssd[[g]] <- setv(sprintf("sigma_slab_%s", g), rprior(dp$slab_sd))
  }
  draw <- function(i) {                                                    # one parameter of the table (shared nodes once)
    r <- pt[i, ]; n <- r$node
    if (is.na(n)) return(r$fixed)
    if (!is.null(truth[[n]])) return(truth[[n]])
    if (!is.na(r$incl)) { g <- if (!is.na(r$group)) r$group else if (r$kind == "loading") "load" else "reg"
      return(setv(n, stats::rbinom(1, 1, pin[[g]]) * stats::rnorm(1, 0, ssd[[g]]))) }
    pr <- r$prior
    if (grepl("^hier\\(", pr)) { f <- sub("^hier\\((.*)\\)$", "\\1", pr); mu <- truth[[sprintf("mu_load_%s", f)]] %||% setv(sprintf("mu_load_%s", f), stats::rnorm(1, 0, sqrt(HIER_V)))
      return(setv(n, exp(stats::rnorm(1, mu, 1)))) }
    pr <- sub("^sd ", "", pr)
    setv(n, rprior(pr))
  }
  val <- function(kind, lhs, op = NULL, rhs = NULL) {
    i <- which(pt$kind == kind & pt$lhs == lhs & (if (is.null(op)) TRUE else pt$op == op) & (if (is.null(rhs)) TRUE else pt$rhs == rhs))
    if (length(i)) draw(i[1]) else 0
  }
  # thresholds and step ratios (drawn as the JAGS code defines them)
  thr <- list(); rat <- list()
  for (v in sp$ordinal) {
    K <- sp$ord$K[match(v, sp$ord$items)]; sv <- safe(v); ty <- sp$ord$type[[v]]
    if (ty == "grm") { b <- rprior(dp$threshold); if (K > 2) for (k in 2:(K - 1)) b <- c(b, b[k - 1] + rprior(dp$threshold_step)) }
    else b <- rprior(dp$threshold, K - 1)
    for (k in seq_len(K - 1)) truth[[sprintf("b_%s[%d]", sv, k)]] <- b[k]
    thr[[v]] <- b; rat[[v]] <- if (ty == "tppcm") c(1, if (K > 2) rprior(dp$step_ratio, K - 2)) else rep(1, K - 1)
  }
  for (i in seq_len(nrow(pt))) if (pt$kind[i] != "threshold" && pt$kind[i] != "step_disc") draw(i)
  for (v in sp$ordinal) if (sp$ord$type[[v]] == "tppcm") {                 # step discriminations: loading x ratio
    K <- sp$ord$K[match(v, sp$ord$items)]; r <- ld[ld$rhs == v & !ld$zero, ]
    for (k in seq_len(K - 2) + 1) truth[[sprintf("as_%s[%d]", safe(v), k)]] <- val("loading", r$lhs, "=~", v) * rat[[v]][k]
  }
  # moderators and person-specific multipliers
  zv <- function(z) as.numeric(data[[z]])
  msum <- function(kind, lhs, rhs_fun) { out <- rep(0, N)
    for (z in sp$modvars) { i <- which(pt$kind == kind & pt$lhs == lhs & pt$rhs == rhs_fun(z)); if (length(i)) out <- out + draw(i[1]) * zv(z) }
    out }
  # persons
  Fv <- list()
  for (f in fs) {
    rr <- pt[pt$kind == "beta" & pt$lhs == f, ]
    mu <- rep(0, N)
    for (k in seq_len(nrow(rr))) if (!rr$rhs[k] %in% fs) mu <- mu + draw(which(pt$kind == "beta" & pt$lhs == f & pt$rhs == rr$rhs[k])[1]) * zv(rr$rhs[k])
    sdf <- if (is.na(sp$fsd[[f]])) val("fsd", f) else sp$fsd[[f]]
    sdj <- sdf * exp(msum("psi", sprintf("V(%s)", f), identity))
    if (f == fs[1] || length(fs) == 1) Fv[[f]] <- mu + sdj * stats::rnorm(N)
    else {
      f1 <- fs[1]
      if (any(rr$rhs == f1)) mu <- mu + draw(which(pt$kind == "beta" & pt$lhs == f & pt$rhs == f1)[1]) * Fv[[f1]]
      if (isTRUE(sp$cov_free)) {
        rho <- val("corr", f1); sd1 <- (if (is.na(sp$fsd[[f1]])) val("fsd", f1) else sp$fsd[[f1]]) * exp(msum("psi", sprintf("V(%s)", f1), identity))
        m1 <- rep(0, N); r1 <- pt[pt$kind == "beta" & pt$lhs == f1, ]
        for (k in seq_len(nrow(r1))) m1 <- m1 + draw(which(pt$kind == "beta" & pt$lhs == f1 & pt$rhs == r1$rhs[k])[1]) * zv(r1$rhs[k])
        Fv[[f]] <- mu + rho * sdj / sd1 * (Fv[[f1]] - m1) + sdj * sqrt(1 - rho^2) * stats::rnorm(N)
      } else Fv[[f]] <- mu + sdj * stats::rnorm(N)
    }
  }
  # indicators
  out <- data
  lam_eff <- function(f, v) val("loading", f, "=~", v) * exp(msum("alpha", sprintf("E(%s)", v), function(z) paste0(z, ":", f)))
  for (v in c(sp$binary, sp$ordinal, sp$cont)) {
    r <- ld[ld$rhs == v & !ld$zero, ]
    dl <- msum("delta", sprintf("E(%s)", v), identity)
    if (v %in% sp$ordinal) {
      a <- lam_eff(r$lhs[1], v); P <- cat_probs(sp$ord$type[[v]], a, Fv[[r$lhs[1]]] - dl, thr[[v]], rat[[v]])
      y <- 1 + rowSums(stats::runif(N) > t(apply(P, 1, cumsum))[, -ncol(P), drop = FALSE])
      y <- sp$ord$levels[[v]][y]
    } else {
      eta <- 0; for (k in seq_len(nrow(r))) eta <- eta + lam_eff(r$lhs[k], v) * Fv[[r$lhs[k]]]
      if (v %in% sp$binary) {
        if (any(dl != 0)) eta <- eta - lam_eff(r$lhs[1], v) * dl
        y <- stats::rbinom(N, 1, stats::plogis(val("intercept", v, "~1") + eta))
      } else {
        sg <- val(pt$kind[pt$lhs == v & pt$op == "~~"][1], v, "~~") * exp(msum("kappa", sprintf("V(%s)", v), identity))
        y <- val("intercept", v, "~1") + eta + dl + sg * stats::rnorm(N)
      }
    }
    y[is.na(data[[v]])] <- NA
    out[[v]] <- y
  }
  list(data = out, truth = truth)
}
