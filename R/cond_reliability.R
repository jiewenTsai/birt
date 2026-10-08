# Conditional reliability of the first latent variable (RTMB engine) ----------------------------
#
# The RTMB engine writes the latent variables as transforms of independent standard normals
# (z1, z2) (see lat_values()). At a value theta of the first latent variable, z1 is fixed and the
# second latent variable (speed) is a nuisance at its conditional median (z2 = 0). The expected
# information of every indicator about (z1, z2) is computed from the parameters that z moves:
#   binary:      eta = d + sum lambda f                      I = P(1 - P) grad(eta) grad(eta)'
#   continuous:  psi = (mu, log sd, skewness)                I = J' I_psi J,  J = d psi / d z
# with I_psi of the normal (diag(1/sd^2, 2)) or of the standardized SHASH (Gauss-Hermite over the
# residual, scores by central differences). J and grad(eta) come from central differences of
# lat_values() and the moderator sums, so cross-loadings, paths, correlations, SHASH latent
# variables and scale / skewness moderation by the latent variables are all included. With the
# prior precision (identity in z), V = [(sum I + I_2)^-1]_11 is the approximate posterior variance
# of z1 and rho(theta) = 1 - V the conditional reliability (Nicewander, 2018: I / (I + 1) when the
# prior variance is 1); speed enters through the Schur complement, so its uncertainty is not
# counted as information about theta.

#' Conditional reliability of the first latent variable
#'
#' For a model fitted with `engine = "rtmb"`, the reliability of the first latent variable (e.g.
#' ability) as a function of its value, \eqn{\rho(\theta) = I(\theta) / (I(\theta) + 1/\sigma^2)}
#' (Nicewander, 2018), where \eqn{I(\theta)} is the expected Fisher information of all indicators
#' about \eqn{\theta}, with the second latent variable (speed) treated as a nuisance (Schur
#' complement, at its conditional median) and \eqn{\sigma^2} the variance of \eqn{\theta} given
#' the covariates. The information is split into the accuracy part (categorical indicators:
#' binary and ordinal), the response-time part (continuous indicators only) and both; for
#' distributional models (`V(y) ~ f` in the syntax, `shash(skew = )` moderated by the latent
#' variables) the scale and the shape of the response-time distribution contribute too.
#'
#' * Ordinal indicators contribute their item information (graded response: Samejima, 1969;
#'   partial credit: the variance of the cumulative step discriminations), the same as
#'   [ordinal_information()].
#' * Moderated loadings and shifts (`E(y) ~ z:f`, `E(y) ~ z`) are evaluated at the moderator
#'   values `at` (default: their means), so the curve is the reliability for persons with these
#'   moderator values.
#' * For a fit with censored indicators (`birt(..., censor = )`) the information is that of the
#'   censored data at the fit's limits: the density part between the limits plus the probability
#'   of a censored value times the squared score of log S(c) or log F(c) (Escobar & Meeker,
#'   1998). [rt_information()] shows what other limits would keep.
#'
#' The curve is checked against the data: `empirical` bins the persons by their EAP and compares
#' \eqn{1 - \mathrm{PSD}^2/\sigma^2} from the adaptive quadrature with the curve at the bin mean.
#' @param object A fit from `birt(..., engine = "rtmb")`.
#' @param theta Values of the first latent variable.
#' @param at Values of the covariates and moderators (a named vector); default their means.
#' @param se Delta-method standard errors of the curves (from `vcov()`).
#' @param bins Number of EAP bins of the empirical check (0: none).
#' @param k Gauss-Hermite nodes over SHASH residuals.
#' @return A data frame (class `birt_cond_rel`) with `theta`, the information `info`, `info_acc`,
#'   `info_rt` and the reliabilities `rel`, `rel_acc`, `rel_rt` (and `se_*` with `se = TRUE`);
#'   attribute `empirical` with the binned check. `plot()` draws the curves.
#' @references Nicewander, W. A. (2018). Conditional reliability coefficients for test scores.
#'   *Psychological Methods, 23*(2), 351-362.
#'
#' Escobar, L. A., & Meeker, W. Q. (1998). Fisher information matrices with censoring,
#' truncation, and explanatory variables. *Statistica Sinica, 8*, 221-237.
#'
#' Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores.
#' *Psychometrika Monograph Supplement, 34*(4, Pt. 2).
#' @examples
#' \donttest{
#' d <- sim_rtirt(N = 500, K = 6)
#' m <- paste("theta =~ y1 + y2 + y3 + y4 + y5 + y6",
#'            "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5 + -1*t6",
#'            "theta =~ t1 + t2 + t3 + t4 + t5 + t6", sep = "\n")
#' fit <- birt(paste(m, "V(t1 + t2 + t3 + t4 + t5 + t6) ~ theta", sep = "\n"), d,
#'             control = rtmb_control(method = "aghq"), progress = FALSE)
#' cr <- cond_reliability(fit)
#' plot(cr)
#' }
#' @export
cond_reliability <- function(object, theta = seq(-3, 3, by = 0.1), at = NULL, se = TRUE, bins = 10, k = 41) {
  if (!inherits(object, "birt_rtmb")) stop("cond_reliability() needs a fit from birt(..., engine = \"rtmb\")", call. = FALSE)
  S <- object$S; r <- object$rtmb
  xb <- cr_at(S, at)
  gh <- statmod::gauss.quad.prob(k, "normal")
  curve <- function(p) cr_curve(p, S, theta, xb, gh)
  out <- curve(r$par)
  if (se) {
    skel <- r$par; x0 <- r$x; V <- r$vcov
    f <- function(x) unlist(curve(vec2list(x, skel))[c("rel", "rel_acc", "rel_rt")])
    G <- vapply(seq_along(x0), function(i) {
      h <- 1e-4 * max(1, abs(x0[i])); e <- replace(numeric(length(x0)), i, h)
      (f(x0 + e) - f(x0 - e)) / (2 * h)
    }, numeric(3 * length(theta)))
    sev <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    n <- length(theta)
    out$se_rel <- sev[seq_len(n)]; out$se_rel_acc <- sev[n + seq_len(n)]; out$se_rel_rt <- sev[2 * n + seq_len(n)]
  }
  emp <- NULL
  if (bins > 0) {
    sc <- object$scores; f1 <- S$fs[1]; s1 <- out$sd_theta[1]
    m <- sc[[f1]]; ps <- sc[[paste0(f1, "_psd")]]
    br <- unique(stats::quantile(m, seq(0, 1, length.out = bins + 1)))
    g <- cut(m, br, include.lowest = TRUE)
    emp <- data.frame(theta = tapply(m, g, mean), rel_empirical = tapply(1 - ps^2 / s1^2, g, mean), n = as.vector(table(g)))
    emp$rel_model <- if (length(theta) > 1) stats::approx(out$theta, out$rel, emp$theta, rule = 2)$y else out$rel
    rownames(emp) <- NULL
  }
  structure(out, class = c("birt_cond_rel", "data.frame"), empirical = emp, factor = S$fs[1], at = xb)
}

# covariate and moderator values: their means, with the values of `at` (a named list or vector)
cr_at <- function(S, at) {
  xb <- if (ncol(S$X)) colMeans(S$X) else numeric(0)
  if (length(at)) {
    bad <- setdiff(names(at), colnames(S$X))
    if (length(bad)) stop("`at`: not covariates of the model: ", paste(bad, collapse = ", "), call. = FALSE)
    xb[names(at)] <- unlist(at)
  }
  xb
}

# the curves at the parameter list p
# the linear predictors of the indicators at theta (speed at its conditional median, covariates
# at xb), their derivatives along z1 and z2, and the information of the binary indicators
cr_setup <- function(p, S, theta, xb) {
  n <- length(theta); S2 <- S
  S2$X <- matrix(xb, n, length(xb), byrow = TRUE, dimnames = list(NULL, colnames(S$X)))
  f1i <- S$f[[1]]
  ls1 <- if (f1i$free_sd) p$lsd1 else log(f1i$sd)
  if (length(f1i$scale)) ls1 <- ls1 + modsum(p$alpha1, f1i$scale, S2$X[1, , drop = FALSE], NULL, S$fs)
  s1 <- exp(ls1)[1]; mu1 <- lin_x(p$beta1, f1i$cols, S2$X[1, , drop = FALSE])[1]
  z1 <- (theta - mu1) / s1; z2 <- rep(0, n); nF <- S$nF
  rows <- seq_len(n); h <- 1e-4
  parts <- function(z1, z2) cr_parts(p, S2, lat_values(p, S2, z1, if (nF == 2) z2 else 0, rows))
  P0 <- parts(z1, z2)
  D1 <- Map(function(a, b) (a - b) / (2 * h), parts(z1 + h, z2), parts(z1 - h, z2))
  D2 <- if (nF == 2) Map(function(a, b) (a - b) / (2 * h), parts(z1, z2 + h), parts(z1, z2 - h)) else lapply(D1, function(v) 0 * v)
  Ib <- cr_zero(n)
  if (length(S$B)) {
    pr <- stats::plogis(P0$eta)
    Ib <- cr_add(Ib, pr * (1 - pr), D1$eta, D2$eta)
  }
  if (length(S$O)) {                       # ordinal: item information about f - delta (cat_info), times d(f - delta)/dz
    Io <- vapply(seq_along(S$O), function(j) cat_info(S$otype[j], P0$A[, j], P0$th[, j], thresholds(p, S, j),
                                                      vapply(seq_len(S$K[j] - 1), function(k) step_ratio(p, S, j, k), 0)), numeric(n))
    Ib <- cr_add(Ib, matrix(Io, n), D1$th, D2$th)
  }
  list(n = n, nF = nF, s1 = s1, P0 = P0, D1 = D1, D2 = D2, Ib = Ib)
}

# information matrices (entries 11, 12, 22 about z1, z2), one value per theta
cr_zero <- function(n) list(i11 = rep(0, n), i12 = rep(0, n), i22 = rep(0, n))
cr_add <- function(acc, w, g1, g2) { acc$i11 <- acc$i11 + rowSums(w * g1 * g1); acc$i12 <- acc$i12 + rowSums(w * g1 * g2); acc$i22 <- acc$i22 + rowSums(w * g2 * g2); acc }

# posterior variance of z1 with the N(0, 1) priors of z1 and z2: reliability and information about theta
cr_rel <- function(I, nF, s1) {
  if (nF == 1) V <- 1 / (1 + I$i11)
  else { a <- 1 + I$i11; b <- I$i12; c <- 1 + I$i22; V <- c / (a * c - b * b) }
  list(rel = 1 - V, info = (1 / V - 1) / s1^2)
}

cr_curve <- function(p, S, theta, xb, gh) {
  st <- cr_setup(p, S, theta, xb); n <- st$n; nF <- st$nF; s1 <- st$s1; P0 <- st$P0; D1 <- st$D1; D2 <- st$D2
  Ib <- st$Ib; Ic <- cr_zero(n); add <- cr_add
  cj <- if (is.null(S$cens)) integer() else S$cens$j
  for (j in seq_along(S$C)) {
    if (j %in% cj) {                                               # censored (censor =): information under the fit's limits
      cs <- S$cens[S$cens$j == j, ]
      Ic <- Map(`+`, Ic, cens_item_info(rt_item(st, p, S, j), st, p, S, j, cs$lower, cs$upper))
      next
    }
    sd <- exp(P0$lsc[, j]); s <- match(j, S$sh)
    g <- list(mu = D1$mu[, j], lsc = D1$lsc[, j], e = D1$e[, j])
    g2 <- list(mu = D2$mu[, j], lsc = D2$lsc[, j], e = D2$e[, j])
    if (is.na(s)) {                                                # normal: I_psi = diag(1/sd^2, 2)
      Ic <- add(Ic, cbind(1 / sd^2, 2), cbind(g$mu, g$lsc), cbind(g2$mu, g2$lsc))
    } else {                                                       # SHASH: scores over the residual
      e <- P0$e[, j]; dl <- exp(p$ldl[s])
      for (q in seq_along(gh$nodes)) {
        sv <- shash_scores(gh$nodes[q], e, dl); x <- sv$x; dx <- sv$dx; de <- sv$de
        s1v <- -dx / sd * g$mu + (-1 - x * dx) * g$lsc + de * g$e       # d log f / d z1
        s2v <- -dx / sd * g2$mu + (-1 - x * dx) * g2$lsc + de * g2$e
        Ic <- add(Ic, matrix(gh$weights[q], n, 1), cbind(s1v), cbind(s2v))
      }
    }
  }
  rel_of <- function(I) cr_rel(I, nF, s1)
  tot <- rel_of(Map(`+`, Ib, Ic)); acc <- rel_of(Ib); rt <- rel_of(Ic)
  data.frame(theta = theta, info = tot$info, info_acc = acc$info, info_rt = rt$info,
             rel = tot$rel, rel_acc = acc$rel, rel_rt = rt$rel, sd_theta = s1)
}

# the standardized SHASH value x at the standard normal quantile z and the derivatives of its log
# density with respect to x and to the skewness (central differences)
shash_scores <- function(z, e, dl, h = 1e-5) {
  x <- T_std(z, e, dl)
  list(x = x, dx = (ld_std(x + h, e, dl) - ld_std(x - h, e, dl)) / (2 * h), de = (ld_std(x, e + h, dl) - ld_std(x, e - h, dl)) / (2 * h))
}

# linear predictors of the indicators at the latent values Fl (rows of S$X): eta (binary), mu,
# log sd and skewness (continuous)
cr_parts <- function(p, S, Fl) {
  X <- S$X; ld <- S$ld; n <- length(Fl[[1]])
  eta <- matrix(0, n, length(S$B)); mu <- lsc <- e <- matrix(0, n, length(S$C)); th <- A <- matrix(0, n, length(S$O))
  # loadings with their moderation (E(y) ~ z:f) and shifts (E(y) ~ z) at the rows of X, as in ind_ll()
  for (j in seq_along(S$B)) {
    v <- rep(p$d[j], n); L1 <- NULL
    for (r in S$ld_b[[j]]) { L <- loading_at(p, S, r, X); if (is.null(L1)) L1 <- L; v <- v + L * Fl[[ld$f[r]]] }
    dl <- mod_eff(p, S, X, S$B[j], "beta"); if (!is.null(dl)) v <- v - L1 * dl
    eta[, j] <- v
  }
  for (j in seq_along(S$O)) {                                       # ordinal: f - delta and the discrimination
    r <- S$ld_o[[j]]; A[, j] <- loading_at(p, S, r, X); v <- Fl[[ld$f[r]]]
    dl <- mod_eff(p, S, X, S$O[j], "beta"); if (!is.null(dl)) v <- v - dl
    th[, j] <- v
  }
  for (j in seq_along(S$C)) {
    v <- rep(p$xi[j], n)
    for (r in S$ld_c[[j]]) v <- v + loading_at(p, S, r, X) * Fl[[ld$f[r]]]
    dl <- mod_eff(p, S, X, S$C[j], "beta"); if (!is.null(dl)) v <- v + dl
    mu[, j] <- v
    kr <- which(S$kap$j == j)
    lsc[, j] <- p$lsig[j] + if (length(kr)) modsum(p$kap[kr], S$kap$m[kr], X, Fl, S$fs) else 0
    s <- match(j, S$sh)
    if (!is.na(s)) {
      er <- which(S$eta$j == j)
      e[, j] <- p$eps[s] + if (length(er)) modsum(p$eta[er], S$eta$m[er], X, Fl, S$fs) else 0
    }
  }
  list(eta = eta, mu = mu, lsc = lsc, e = e, th = th, A = A)
}

#' @rdname cond_reliability
#' @param x A `birt_cond_rel` object.
#' @param ... Passed to [graphics::plot()].
#' @export
plot.birt_cond_rel <- function(x, ...) {
  f <- attr(x, "factor")
  graphics::plot(x$theta, x$rel, type = "l", lwd = 2, ylim = c(0, 1), xlab = f, ylab = "conditional reliability", ...)
  graphics::lines(x$theta, x$rel_acc, lty = 2, lwd = 2, col = "grey40")
  graphics::lines(x$theta, x$rel_rt, lty = 3, lwd = 2, col = "grey40")
  if (!is.null(x$se_rel)) {
    graphics::lines(x$theta, pmin(x$rel + 1.96 * x$se_rel, 1), lty = 1, col = "grey70")
    graphics::lines(x$theta, pmax(x$rel - 1.96 * x$se_rel, 0), lty = 1, col = "grey70")
  }
  emp <- attr(x, "empirical")
  if (!is.null(emp)) graphics::points(emp$theta, emp$rel_empirical, pch = 19)
  graphics::legend("bottom", c("accuracy + RT", "accuracy only", "RT only", if (!is.null(emp)) "EAP bins (AGHQ)"),
                   lty = c(1, 2, 3, NA), pch = c(NA, NA, NA, if (!is.null(emp)) 19), col = c("black", "grey40", "grey40", "black"),
                   lwd = 2, bty = "n", horiz = TRUE, cex = 0.8)
  invisible(x)
}

#' @export
print.birt_cond_rel <- function(x, digits = 3, ...) {
  cat(sprintf("Conditional reliability of %s (accuracy + RT, accuracy only, RT only)\n", attr(x, "factor")))
  y <- as.data.frame(unclass(x))[, intersect(c("theta", "rel", "se_rel", "rel_acc", "rel_rt", "info"), names(x))]
  i <- unique(round(seq(1, nrow(y), length.out = min(nrow(y), 13))))
  print(format(y[i, ], digits = digits), row.names = FALSE)
  emp <- attr(x, "empirical")
  if (!is.null(emp)) { cat("\nCheck: persons binned by EAP, 1 - PSD^2 / var (AGHQ) against the curve\n"); print(format(emp, digits = digits), row.names = FALSE) }
  invisible(x)
}
