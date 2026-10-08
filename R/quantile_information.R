# Information of each quantile regression of the response times about the first latent variable ----
#
# The tau-th regression of item i, q_i(tau | speed, theta), is estimated by the check loss, whose
# estimating function is psi = tau - 1{t <= q}. Its information about z = (z1, z2) (Godambe, 1960)
# is E[d psi / dz] E[d psi / dz]' / Var(psi) = f(q)^2 (dq/dz)(dq/dz)' / (tau (1 - tau)), the
# information of the indicator 1{t <= q}: with g = d log F(c) / dz at c = q (rt_item()$tail_grad),
# tau / (1 - tau) g g'. Rank one per item; summed over the items, speed is projected out with its
# N(0, 1) prior as in cond_reliability() (cr_rel()).

#' Information of each quantile regression of the response times
#'
#' For a model fitted with `engine = "rtmb"`: how much the \eqn{\tau}-th quantile regression of the
#' response times on the latent variables (e.g. the cross relation `t ~ speed + ability`, with
#' cross-loadings `ability =~ t1` and scale or skew moderation `V(t1) ~ ability`) says about the
#' first latent variable, and which \eqn{\tau} says the most. Speed is a nuisance and is projected
#' out, as in [cond_reliability()].
#'
#' The \eqn{\tau}-th regression of item i is estimated by the check loss, whose estimating function
#' is \eqn{\psi_\tau = \tau - 1\{t_{ij} \le q_{i\tau}\}}. Its information about the latent variable
#' (Godambe, 1960) is
#' \deqn{I_{i\tau}(\theta) = \frac{f_i(q_{i\tau})^2\, \beta_i(\tau)^2}{\tau (1 - \tau)},
#'   \qquad \beta_i(\tau) = \partial q_{i\tau} / \partial\theta,}
#' a signal (how far the quantile moves with \eqn{\theta}) times the precision of the quantile, the
#' inverse of its asymptotic variance \eqn{\tau(1-\tau)/f(q_\tau)^2} (Koenker, 2005). It equals the
#' Fisher information of the indicator "faster than \eqn{q_{i\tau}}". The items are summed (with
#' the speed direction, item by item rank one) and speed is projected out.
#'
#' The values for different \eqn{\tau} are not additive: each is the information of one regression
#' on its own, and at most the information of the whole response-time distribution (`info_rt`, as in
#' [cond_reliability()]); `efficiency` is the ratio. For one item with a normal residual and a
#' location effect (a cross-loading) the median regression is best and keeps \eqn{2/\pi} = .64 of
#' the item's information; for a scale effect (`V(t) ~ ability`) the median carries nothing and the
#' best single quantiles are \eqn{\tau} = .058 and .942, which keep .30 (maximizing
#' \eqn{z_\tau^2 \phi(z_\tau)^2 / (\tau(1-\tau))}; on the choice of sample quantiles, Ogawa, 1951).
#' With several items and the projection of speed, the efficiencies of the curve differ somewhat
#' from these item values.
#'
#' The quantities are evaluated at the parameter estimates (ML or ELGM), at each `theta` with speed at
#' its conditional median and the covariates at `at`; the effects are conditional on \eqn{\theta}, not
#' the average marginal effects of [quantiles()]. The response-time distribution is that of the
#' fitted model without censoring.
#' @param object A fit of `birt(..., engine = "rtmb")` with continuous indicators.
#' @param p Quantile levels \eqn{\tau}.
#' @param theta Values of the first latent variable.
#' @param at Covariate and moderator values (default: their means), as in [cond_reliability()].
#' @return A list of class `birt_qinfo`:
#'   * `curve`: `theta`, `p`, `info` (all items, speed projected out), `info_rt` (the whole RT
#'     distribution), `efficiency`;
#'   * `items`: `theta`, `item`, `p`, `q` (the quantile, model scale), `beta`
#'     (\eqn{\partial q / \partial\theta}), `density` (\eqn{f(q)}), `info` (the item's
#'     \eqn{I_{i\tau}}, speed held fixed);
#'   * `best`: for each `theta`, the \eqn{\tau} of `p` with the largest `info`.
#' @references
#' Godambe, V. P. (1960). An optimum property of regular maximum likelihood estimation.
#' *The Annals of Mathematical Statistics, 31*, 1208-1211. \doi{10.1214/aoms/1177705693}
#'
#' Koenker, R. (2005). *Quantile regression*. Cambridge University Press.
#' \doi{10.1017/CBO9780511754098}
#'
#' Ogawa, J. (1951). Contributions to the theory of systematic statistics, I. *Osaka Mathematical
#' Journal, 3*, 175-213.
#' @examples
#' \donttest{
#' d <- sim_rtirt(N = 500, K = 6)
#' m <- paste("theta =~ y1 + y2 + y3 + y4 + y5 + y6",
#'            "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5 + -1*t6",
#'            "V(t1 + t2 + t3 + t4 + t5 + t6) ~ theta", sep = "\n")
#' fit <- birt(m, d, progress = FALSE)
#' qi <- quantile_information(fit)
#' qi
#' plot(qi)
#' }
#' @export
quantile_information <- function(object, p = seq(0.05, 0.95, by = 0.05), theta = c(-2, 0, 2), at = NULL) {
  if (!inherits(object, "birt_rtmb")) stop("quantile_information() needs a fit from birt(..., engine = \"rtmb\")", call. = FALSE)
  S <- object$S; par <- object$rtmb$par
  if (!length(S$C)) stop("the model has no continuous indicators (response times)", call. = FALSE)
  if (!is.numeric(p) || !length(p) || any(p <= 0 | p >= 1)) stop("p must be quantile levels in (0, 1)", call. = FALSE)
  p <- round(p, 12)                                                  # seq() artefacts, so that p %in% c(.75, .9) matches
  xb <- cr_at(S, at)
  st <- cr_setup(par, S, theta, xb); n <- st$n; s1 <- st$s1
  item <- lapply(seq_along(S$C), function(j) rt_item(st, par, S, j))
  info_rt <- rt_information(object, theta, bands = c(0, 1), limits = list(), at = at)$check$info_rt
  curve <- items <- list()
  for (u in p) {
    z <- stats::qnorm(u); I <- cr_zero(n)
    for (j in seq_along(S$C)) {
      g <- item[[j]]$tail_grad(z, upper = FALSE)                       # d log F(q) / d (z1, z2)
      I <- cr_add(I, matrix(u / (1 - u), n, 1), g[, 1, drop = FALSE], g[, 2, drop = FALSE])
      sd <- exp(st$P0$lsc[, j]); s <- match(j, S$sh)
      if (is.na(s)) { x <- rep(z, n); dens <- stats::dnorm(z) / sd }
      else { e <- st$P0$e[, j]; dl <- exp(par$ldl[s]); x <- T_std(rep(z, n), e, dl); dens <- exp(ld_std(x, e, dl)) / sd }
      items[[length(items) + 1]] <- data.frame(theta = theta, item = S$C[j], p = u, q = st$P0$mu[, j] + sd * x,
                                               beta = -g[, 1] * u / dens / s1, density = dens,
                                               info = u / (1 - u) * g[, 1]^2 / s1^2)
    }
    info <- cr_rel(I, st$nF, s1)$info
    curve[[length(curve) + 1]] <- data.frame(theta = theta, p = u, info = info, info_rt = info_rt, efficiency = info / info_rt)
  }
  curve <- do.call(rbind, curve); items <- do.call(rbind, items); rownames(items) <- NULL
  best <- do.call(rbind, lapply(split(curve, curve$theta), function(d) d[which.max(d$info), c("theta", "p", "info", "efficiency")]))
  rownames(best) <- NULL
  structure(list(curve = curve, items = items, best = best), class = "birt_qinfo", factor = S$fs[1])
}

#' @rdname quantile_information
#' @param x A `birt_qinfo` object.
#' @param ... Unused.
#' @export
print.birt_qinfo <- function(x, ...) {
  f <- attr(x, "factor")
  cat(sprintf("Information of the tau-th quantile regression of the response times about %s (speed projected out)\n\n", f))
  w <- stats::reshape(x$curve[, c("theta", "p", "efficiency")], idvar = "p", timevar = "theta", direction = "wide")
  names(w) <- c("p", sprintf("eff at %s = %s", f, format(unique(x$curve$theta))))
  print_table(w, 3)
  cat(sprintf("\nefficiency: information of the regression at p / information of the whole RT distribution\nBest p:\n"))
  print_table(x$best, 3)
  invisible(x)
}

#' @rdname quantile_information
#' @export
plot.birt_qinfo <- function(x, ...) {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1)); on.exit(graphics::par(op))
  d <- x$curve; th <- unique(d$theta); f <- attr(x, "factor")
  panel <- function(y, ylim, ylab, main) {
    graphics::plot(NA, xlim = c(0, 1), ylim = ylim, xlab = expression("quantile level " * tau), ylab = ylab, main = main)
    for (i in seq_along(th)) graphics::lines(d$p[d$theta == th[i]], y[d$theta == th[i]], type = "b", pch = 16, lwd = 2, col = i)
  }
  panel(d$info, c(0, max(d$info)), sprintf("information about %s", f), "Each quantile regression")
  graphics::legend("top", sprintf("%s = %s", f, format(th)), col = seq_along(th), lwd = 2, bty = "n")
  panel(d$efficiency, c(0, 1), "efficiency (share of the whole RT information)", "Against the whole RT distribution")
  invisible(x)
}
