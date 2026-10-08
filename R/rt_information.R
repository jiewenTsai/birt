# Where in the response-time distribution the information about the first latent variable lies --------
#
# At a value theta of the first latent variable (speed a nuisance at its conditional median, as in
# cond_reliability()), the score of a log response time t_i about (z1, z2) at its quantile level u is
# s(u) = d log f(Q_i(u)) / d z. Integrating over u in (0, 1) gives the information of cond_reliability().
# * partition: the efficient score about z1, s1 - c s2 with c = I12 / (1 + I22) (speed projected out
#   with its N(0, 1) prior), squared and summed over the items, is a density over u that integrates to
#   the RT information about theta (up to the small prior term c^2 of speed);
# * designs: with a time limit at the item quantile u_c the observations above it are right censored
#   (only "slower than c" is known), so the information is int_0^{u_c} s s' du + (1 - u_c) g g',
#   g = d log S(c) / d z; with a limit u_f for fast responses (left censored) it is
#   u_f g g' + int_{u_f}^1 s s' du, g = d log F(c) / d z (Efron & Johnstone, 1990; Escobar & Meeker, 1998).

#' Where response times carry information about the first latent variable
#'
#' For a model fitted with `engine = "rtmb"`: how the information that the response times carry about
#' the first latent variable (e.g. ability, through cross-loadings `ability =~ t1` and scale or skew
#' moderation `V(t1) ~ ability`) is spread over the quantile levels of the response-time
#' distribution, and how much of it, and of the reliability, a time limit keeps. Speed is a nuisance
#' and is projected out, as in [cond_reliability()].
#'
#' * `bands`: the information of the response times in each band of quantile levels (all items at
#'   their own quantiles). The bands partition the RT information of [cond_reliability()]
#'   (`info_rt`); they show where it lies, not what a design would keep.
#' * `limits`: the information and reliability when every response time beyond a limit is censored
#'   (only known to exceed it): `side = "slow"`, a time limit at the item's quantile u (as with a
#'   time limit); `side = "fast"`, responses faster than the quantile u known only as fast (as when
#'   rapid responses are not modelled). Censoring keeps the information that a response was beyond
#'   the limit, so a slow tail can hold much of the partition yet lose little under a time limit.
#'   Discarding the responses instead (truncation) loses more and biases a model that ignores it.
#'
#' The decomposition of Fisher information over an ordered sample follows Zheng & Gastwirth (2000);
#' censored information Efron & Johnstone (1990) and Escobar & Meeker (1998). The curves describe
#' the fitted model: fast responses that are rapid guesses (Wise & Kong, 2005) violate it, and their
#' apparent information is then an artefact.
#' The response-time distribution is that of the fitted model without censoring; for a fit with
#' `censor =`, [cond_reliability()] gives the information under the fit's own limits. Binary and
#' ordinal indicators enter `rel` as in [cond_reliability()], and moderated loadings or shifts
#' are evaluated at the moderator values `at`.
#' @param object A fit of `birt(..., engine = "rtmb")` with continuous indicators.
#' @param theta Values of the first latent variable.
#' @param bands Breaks of the quantile bands.
#' @param limits Quantile levels of the limits; `slow` for time limits, `fast` for fast responses.
#' @param at Covariate and moderator values (default: their means), as in [cond_reliability()].
#' @param k Gauss-Legendre nodes per interval.
#' @return A list of class `birt_rt_info`: `density` (`theta`, `u`, `info`: the information density over
#'   u), `bands` (`theta`, `band`, `info`, `share`), `limits` (`theta`, `side`, `u`, `info_rt`, `kept`
#'   (share of the RT information), `rel` (reliability with the accuracy items), `rel_full`), and
#'   `check` (the RT information of the partition and of [cond_reliability()]).
#' @references
#' Zheng, G., & Gastwirth, J. L. (2000). Where is the Fisher information in an ordered sample?
#' *Statistica Sinica, 10*, 1267-1280.
#'
#' Efron, B., & Johnstone, I. M. (1990). Fisher's information in terms of the hazard rate.
#' *The Annals of Statistics, 18*, 38-62. \doi{10.1214/aos/1176347492}
#'
#' Escobar, L. A., & Meeker, W. Q. (1998). Fisher information matrices with censoring, truncation,
#' and explanatory variables. *Statistica Sinica, 8*, 221-237.
#'
#' Wise, S. L., & Kong, X. (2005). Response time effort: A new measure of examinee motivation in
#' computer-based tests. *Applied Measurement in Education, 18*, 163-183. \doi{10.1207/s15324818ame1802_2}
#' @export
rt_information <- function(object, theta = c(-2, 0, 2), bands = c(0, .05, .1, .25, .5, .75, .9, .95, 1),
                           limits = list(slow = c(.75, .9, .95, .99), fast = c(.01, .05, .1)), at = NULL, k = 60) {
  if (!inherits(object, "birt_rtmb")) stop("rt_information() needs a fit from birt(..., engine = \"rtmb\")", call. = FALSE)
  S <- object$S; p <- object$rtmb$par
  if (!length(S$C)) stop("the model has no continuous indicators (response times)", call. = FALSE)
  xb <- cr_at(S, at)
  st <- cr_setup(p, S, theta, xb); n <- st$n; s1 <- st$s1
  gl <- statmod::gauss.quad(k, "legendre")
  zlim <- 10
  item <- lapply(seq_along(S$C), function(j) rt_item(st, p, S, j))
  # int over z in (a, b) of phi(z) * s s' summed over the items: entries 11, 12, 22 (n each)
  Iz <- function(a, b) {
    a <- max(a, -zlim); b <- min(b, zlim); out <- cr_zero(n)
    if (b <= a) return(out)
    zq <- (b - a) / 2 * gl$nodes + (a + b) / 2; wq <- (b - a) / 2 * gl$weights * stats::dnorm(zq)
    for (it in item) for (q in seq_along(zq)) { M <- it$score(zq[q]); out <- cr_add(out, matrix(wq[q], n, 1), M[, 1, drop = FALSE], M[, 2, drop = FALSE]) }
    out
  }
  plus <- function(A, B) Map(`+`, A, B)
  qz <- function(u) stats::qnorm(pmin(pmax(u, 0), 1))
  # full RT information (pieces between the band breaks, so that the partition uses the same nodes)
  br <- sort(unique(c(0, bands, 1)))
  pieces <- lapply(seq_len(length(br) - 1), function(b) Iz(qz(br[b]), qz(br[b + 1])))
  Ic <- Reduce(plus, pieces)
  cc <- Ic$i12 / (1 + Ic$i22)                                         # projection of speed (with its prior)
  eff <- function(I) I$i11 - 2 * cc * I$i12 + cc^2 * I$i22             # efficient part of a piece
  rt_full <- cr_rel(Ic, st$nF, s1); tot_full <- cr_rel(plus(st$Ib, Ic), st$nF, s1)
  band_info <- vapply(pieces, eff, numeric(n)); if (is.null(dim(band_info))) band_info <- matrix(band_info, n)
  lab <- sprintf("(%s, %s]", format(br[-length(br)]), format(br[-1]))
  bands_df <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(theta = theta[i], band = lab, info = band_info[i, ] / s1^2,
                                                                        share = band_info[i, ] / sum(band_info[i, ]))))
  check <- data.frame(theta = theta, info_partition = rowSums(band_info) / s1^2, speed_prior = cc^2 / s1^2, info_rt = rt_full$info)
  # censoring designs
  lim <- list()
  for (side in intersect(names(limits), c("slow", "fast"))) for (u in limits[[side]]) {
    zc <- stats::qnorm(u)
    obs <- if (side == "slow") Iz(-Inf, zc) else Iz(zc, Inf)
    tl <- cr_zero(n)
    for (it in item) { g <- it$tail_grad(zc, upper = side == "slow"); tl <- cr_add(tl, matrix(if (side == "slow") 1 - u else u, n, 1), g[, 1, drop = FALSE], g[, 2, drop = FALSE]) }
    Iu <- plus(obs, tl); r1 <- cr_rel(Iu, st$nF, s1); r2 <- cr_rel(plus(st$Ib, Iu), st$nF, s1)
    lim[[length(lim) + 1]] <- data.frame(theta = theta, side = side, u = u, info_rt = r1$info, kept = r1$info / rt_full$info,
                                         rel = r2$rel, rel_full = tot_full$rel)
  }
  # density over u for plots
  ug <- c(seq(0.002, 0.05, by = 0.002), seq(0.06, 0.94, by = 0.01), seq(0.95, 0.998, by = 0.002))
  dens <- vapply(ug, function(u) { M <- lapply(item, function(it) it$score(stats::qnorm(u)))
    Reduce(`+`, lapply(M, function(m) (m[, 1] - cc * m[, 2])^2)) / s1^2 }, numeric(n))
  if (is.null(dim(dens))) dens <- matrix(dens, n)
  density <- data.frame(theta = rep(theta, length(ug)), u = rep(ug, each = n), info = as.vector(dens))
  structure(list(density = density, bands = bands_df, limits = if (length(lim)) do.call(rbind, lim), check = check),
            class = "birt_rt_info", factor = S$fs[1])
}

#' @rdname rt_information
#' @param x A `birt_rt_info` object.
#' @param ... Unused.
#' @export
print.birt_rt_info <- function(x, ...) {
  f <- attr(x, "factor")
  cat(sprintf("Information of the response times about %s by quantile level of the RT (speed projected out)\n\n", f))
  b <- x$bands; b$share <- round(100 * b$share, 1)
  w <- stats::reshape(b[, c("theta", "band", "share")], idvar = "band", timevar = "theta", direction = "wide")
  names(w) <- sub("^share\\.", sprintf("%% at %s = ", f), names(w))
  print(w, row.names = FALSE)
  if (!is.null(x$limits)) {
    cat("\nCensoring beyond a limit (slow: time limit at the item's quantile u; fast: faster than u known only as fast)\n")
    l <- x$limits; l$kept <- round(l$kept, 3); l$rel <- round(l$rel, 3); l$rel_full <- round(l$rel_full, 3); l$info_rt <- round(l$info_rt, 3)
    print(l, row.names = FALSE)
  }
  invisible(x)
}

#' @rdname rt_information
#' @export
plot.birt_rt_info <- function(x, ...) {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1)); on.exit(graphics::par(op))
  d <- x$density; th <- unique(d$theta); f <- attr(x, "factor")
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, max(d$info)), xlab = "quantile level u of the response time",
                 ylab = sprintf("information about %s per unit of u", f), main = "Where the RT information lies")
  for (i in seq_along(th)) graphics::lines(d$u[d$theta == th[i]], d$info[d$theta == th[i]], lwd = 2, col = i)
  graphics::legend("top", sprintf("%s = %s", f, format(th)), col = seq_along(th), lwd = 2, bty = "n")
  l <- x$limits
  if (!is.null(l)) {
    s <- l[l$side == "slow", ]
    graphics::plot(NA, xlim = range(s$u), ylim = c(min(s$kept, 0.5), 1), xlab = "time limit at the item quantile u",
                   ylab = "share of the RT information kept", main = "Time limits (right censoring)")
    for (i in seq_along(th)) graphics::lines(s$u[s$theta == th[i]], s$kept[s$theta == th[i]], type = "b", lwd = 2, col = i, pch = 16)
  }
  invisible(x)
}

# scores of continuous indicator j about z1 and z2 at the standard normal quantile z of its residual
# (z: one value, or one per theta; n x 2), and d log S(c) / d z, d log F(c) / d z at the time c of
# the quantile z (st: cr_setup())
rt_item <- function(st, p, S, j) {
  n <- st$n
  sd <- exp(st$P0$lsc[, j]); s <- match(j, S$sh); e <- st$P0$e[, j]; dl <- if (is.na(s)) 1 else exp(p$ldl[s])
  g1 <- list(mu = st$D1$mu[, j], lsc = st$D1$lsc[, j], e = st$D1$e[, j]); g2 <- list(mu = st$D2$mu[, j], lsc = st$D2$lsc[, j], e = st$D2$e[, j])
  score <- function(z) {
    if (is.na(s)) { x <- rep_len(z, n); dx <- -x; de <- 0 }
    else { sv <- shash_scores(rep_len(z, n), e, dl); x <- sv$x; dx <- sv$dx; de <- sv$de }
    sc <- function(g) -dx / sd * g$mu + (-1 - x * dx) * g$lsc + de * g$e
    cbind(sc(g1), sc(g2))
  }
  tail_grad <- function(z, upper) {
    mu0 <- st$P0$mu[, j]; ls0 <- st$P0$lsc[, j]; c0 <- mu0 + exp(ls0) * (if (is.na(s)) rep_len(z, n) else T_std(rep_len(z, n), e, dl))
    lp <- function(mu, ls, ee) {
      y <- (c0 - mu) / exp(ls)
      u <- if (is.na(s)) y else { a <- eff_a(ee, dl); sinh(dl * (asinh(y * cosh(a) / dl + sinh(a)) - a)) }
      stats::pnorm(u, lower.tail = !upper, log.p = TRUE)
    }
    h <- 1e-5
    d <- function(g) (lp(mu0 + h * g$mu, ls0 + h * g$lsc, e + h * g$e) - lp(mu0 - h * g$mu, ls0 - h * g$lsc, e - h * g$e)) / (2 * h)
    cbind(d(g1), d(g2))
  }
  # standard normal quantile of the time c (data scale) for each theta
  zof <- function(c) { x <- (c - st$P0$mu[, j]) / sd; if (is.na(s)) x else T_inv(x, e, dl) }
  list(score = score, tail_grad = tail_grad, zof = zof)
}

# information (entries 11, 12, 22 about z1, z2) of continuous indicator j censored below `lower`
# and above `upper` (data scale; -Inf / Inf: none): the density part between the limits plus
# F(l) g_l g_l' + S(c) g_c g_c' (Escobar & Meeker, 1998), with per-theta limits on the z scale
cens_item_info <- function(it, st, p, S, j, lower, upper, k = 60, zlim = 10) {
  n <- st$n; gl <- statmod::gauss.quad(k, "legendre"); out <- cr_zero(n)
  zl <- if (is.finite(lower)) it$zof(lower) else rep(-Inf, n); zu <- if (is.finite(upper)) it$zof(upper) else rep(Inf, n)
  a <- pmax(zl, -zlim); b <- pmax(pmin(zu, zlim), a)
  for (q in seq_along(gl$nodes)) {
    zq <- (b - a) / 2 * gl$nodes[q] + (a + b) / 2; wq <- (b - a) / 2 * gl$weights[q] * stats::dnorm(zq)
    M <- it$score(zq); out <- cr_add(out, matrix(wq, n, 1), M[, 1, drop = FALSE], M[, 2, drop = FALSE])
  }
  if (is.finite(lower)) { g <- it$tail_grad(zl, upper = FALSE); out <- cr_add(out, matrix(stats::pnorm(zl), n, 1), g[, 1, drop = FALSE], g[, 2, drop = FALSE]) }
  if (is.finite(upper)) { g <- it$tail_grad(zu, upper = TRUE); out <- cr_add(out, matrix(stats::pnorm(zu, lower.tail = FALSE), n, 1), g[, 1, drop = FALSE], g[, 2, drop = FALSE]) }
  out
}
