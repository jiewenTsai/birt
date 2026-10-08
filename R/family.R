#' Residual distributions for latent variables and continuous indicators
#'
#' Used in the `family` argument of [birt()] to replace the default normal
#' distribution of a latent variable (its residual, in a latent regression) or of a
#' continuous indicator by an asymmetric Laplace distribution (ALD), which turns the
#' corresponding regression into a quantile regression (Yu & Moyeed, 2001).
#'
#' @param p Quantile level(s) in (0, 1). A vector (for one entry of `family` only)
#'   fits one model per quantile and returns a `birt_quantile` object.
#' @return An object of class `birt_family`.
#' @details
#' `family = list(speed = ald(0.25))` makes the 0.25 quantile of `speed` given its
#' predictors equal to the linear predictor (a latent quantile regression; Wang, Feng &
#' Song, 2016; Zhu, Gao & Zhang, 2021). `family = list(.continuous = ald(0.5))` gives
#' every continuous indicator a median (least absolute deviation) measurement model,
#' which is robust to outlying values such as rapid guesses. Names are resolved by their
#' role in the model syntax: a latent variable name is the structural level, an
#' indicator name the measurement level, and `.continuous` all
#' continuous indicators. Binary indicators always use the logit link.
#'
#' ALD on a latent variable whose variance is fixed (e.g. `theta`, variance 1 for
#' identification) fixes the ALD scale so that the implied variance equals the fixed
#' value; with a free variance the scale is estimated. The factor correlation `f1 ~~ f2`
#' cannot be combined with ALD; write the relation as a regression `f2 ~ f1`.
#'
#' Sampling uses the normal-exponential mixture of Kozumi and Kobayashi (2011). Posterior
#' intervals under a working ALD likelihood can be too narrow (Yang, Wang & He, 2016).
#' @examples
#' ald(0.25)
#' ald(c(0.1, 0.5, 0.9))
#' @export
ald <- function(p = 0.5) {
  if (!is.numeric(p) || !length(p) || any(is.na(p)) || any(p <= 0 | p >= 1))
    stop("ald(): p must be in (0, 1)")
  structure(list(family = "ald", p = p), class = "birt_family")
}

#' @rdname ald
#' @details `normal()` is the default; a heteroscedastic residual is written in the model
#'   syntax, `V(t1) ~ theta` (log residual SD linear in moderators).
#' @export
normal <- function() {
  structure(list(family = "normal", p = NA_real_, scale = NULL), class = "birt_family")
}

#' Sinh-arcsinh (SHASH) residual distribution
#'
#' For `birt(..., engine = "rtmb")`: replaces the normal residual of the second latent
#' variable (e.g. `speed` in `speed ~ theta + x`) or of continuous indicators (e.g. log
#' response times, `family = list(.continuous = shash(skew = "theta"))`) by the
#' sinh-arcsinh distribution of Jones and Pewsey (2009): `e = sinh((asinh(z) + eps) / delta)`
#' with `z ~ N(0, 1)`, where `eps` sets the skewness (`eps > 0`: right-skewed) and `delta` the
#' tail weight (`delta < 1`: heavier tails than the normal). Moderators make the skewness
#' (`skew`) linear in latent variables or covariates, and `V(y) ~ z` in the model syntax the
#' log scale, so that the
#' predictor effects on the conditional quantiles can change with the level p; see
#' [quantiles()]. This is the likelihood-based counterpart of fitting [ald()] at several
#' quantiles with `engine = "jags"`: one fit gives every quantile, and the quantile curves
#' cannot cross.
#' @param skew Moderators of the skewness `eps`.
#' @return An object of class `birt_family`.
#' @examples
#' shash(skew = "theta")
#' @export
shash <- function(skew = NULL) {
  check_mods(skew, "shash(skew = )")
  structure(list(family = "shash", p = NA_real_, scale = NULL, skew = skew), class = "birt_family")
}

check_mods <- function(m, what) if (!is.null(m) && (!is.character(m) || anyNA(m) || any(m == "")))
  stop(what, " takes names of latent variables or covariates", call. = FALSE)

#' @export
print.birt_family <- function(x, ...) {
  cat(fam_label(x), "\n")
  invisible(x)
}

fam_label <- function(f) {
  m <- function(v, what) if (length(v)) sprintf("%s ~ %s", what, paste(v, collapse = " + "))
  extra <- c(m(f$scale, "log sd"), m(f$skew, "skew"))
  base <- switch(f$family, normal = "normal", shash = "SHASH", ald = sprintf("ALD(p = %s)", format(f$p)))
  if (length(extra)) sprintf("%s (%s)", base, paste(extra, collapse = "; ")) else base
}
is_q <- function(f) f$family == "ald"
needs_rtmb <- function(f, fs = character()) f$family == "shash" || any(f$scale %in% fs) || length(f$skew) > 0   # latent scale moderators: RTMB

# Resolve `family` against the model: one entry per factor and per continuous indicator.
resolve_family <- function(family, spec) {
  fac <- stats::setNames(rep(list(normal()), length(spec$factors)), spec$factors)
  ind <- stats::setNames(rep(list(normal()), length(spec$cont)), spec$cont)
  if (is.null(family)) return(list(factor = fac, ind = ind))
  if (!is.list(family) || inherits(family, "birt_family") || is.null(names(family)) || any(names(family) == ""))
    stop("`family` must be a named list, e.g. list(speed = ald(0.25))")
  for (n in names(family)) {
    f <- family[[n]]
    if (!inherits(f, "birt_family")) stop("family$", n, " must be created by ald(), normal() or shash()")
    if (length(f$p) > 1) stop("internal: one quantile per fit")
  }
  all_cont <- intersect(".continuous", names(family))
  if (length(all_cont)) for (v in spec$cont) ind[[v]] <- family[[all_cont]]
  for (n in setdiff(names(family), all_cont)) {
    if (n %in% spec$factors) fac[[n]] <- family[[n]]
    else if (n %in% spec$cont) ind[[n]] <- family[[n]]
    else if (n %in% (spec$ordinal %||% character()))
      stop("family$", n, ": the model of an ordinal indicator is set by itemtype =", call. = FALSE)
    else if (n %in% spec$binary)
      stop("family$", n, ": binary indicators use the logit link; ald() applies to continuous indicators and latent variables")
    else stop("family: '", n, "' is not in the model. Latent variables: ",
              paste(spec$factors, collapse = ", "), "; continuous indicators: ",
              if (length(spec$cont)) paste(spec$cont, collapse = ", ") else "(none)")
  }
  ald_f <- names(fac)[vapply(fac, is_q, logical(1))]
  if (length(ald_f) && spec$cov_free)
    stop("ald() on ", paste(ald_f, collapse = ", "), " cannot be combined with the factor correlation '",
         spec$factors[1], " ~~ ", spec$factors[2], "'; write it as a regression, e.g. '",
         spec$factors[2], " ~ ", spec$factors[1], "'")
  list(factor = fac, ind = ind)
}
