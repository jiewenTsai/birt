# Influence functions (casewise scores) of quantile effects -----------------------------------

#' Influence functions of quantile effects
#'
#' For a maximum likelihood fit of [birt()] with `engine = "rtmb"`: the casewise influence
#' function of each quantile effect of [quantiles()] (e.g. the cross-relation `t ~ speed + ability`
#' of a log response time at p = 0.9). A quantile effect is a smooth function
#' \eqn{\beta(p) = g_p(\vartheta)} of the parameters, so to first order
#' \deqn{\hat\beta(p) - \beta(p) \approx \frac1N \sum_j \mathrm{IF}_j, \qquad
#'   \mathrm{IF}_j = N\, \nabla g_p^\top \hat V s_j,}
#' with \eqn{s_j} the person's score ([estfun.birt_rtmb()]) and \eqn{\hat V} the inverse observed
#' information: the influence function of a function of an M-estimator (Hampel, 1974), obtained by
#' the delta method. It gives
#' * robust standard errors, \eqn{\{\sum_j \mathrm{IF}_j^2\}^{1/2} / N} (the sandwich of Huber,
#'   1967, and White, 1982, valid when the SHASH or normal distribution is misspecified), next to
#'   the model-based ones of [quantiles()];
#' * the influence of each person on an effect: removing person j changes the estimate by about
#'   \eqn{-\mathrm{IF}_j / N} (as in Cook, 1977);
#' * score-based tests of the invariance of an effect along a moderator, `score_test(qs, z)`
#'   (Zeileis & Hornik, 2007; Merkle & Zeileis, 2013).
#'
#' **Conditional or unconditional.** The effects are averages over the persons' covariates
#' (average marginal effects). `type = "conditional"` treats the covariates as fixed, as
#' [quantiles()] does. `type = "unconditional"` also counts their sampling: \eqn{\mathrm{IF}_j}
#' gains \eqn{h_j - \hat\beta(p)}, the deviation of person j's own marginal effect from the
#' average (Graubard & Korn, 1999). The two differ only when the effect depends on covariates.
#'
#' **ELGM fits.** `quantile_score()` is a maximum likelihood tool and refuses an ELGM fit (the
#' default of [birt()]): refit the same model with `control = rtmb_control(method = "aghq")`.
#' Influence functions and score tests need casewise scores that sum to zero at the estimate; at
#' the ELGM posterior mode they sum to minus the gradient of the log prior, the reported ELGM
#' estimates are posterior means rather than the mode, and MAP versions of the score tests need
#' centred scores, with error rates that depend on the covariance estimate (Debelak, Pawel,
#' Strobl & Merkle, 2022). The Bayesian counterparts: [quantiles()] of the ELGM fit gives
#' posterior means, SDs and intervals of the quantile effects, [quantile_test()] works on ELGM
#' fits, and the invariance of an effect along a numeric `z` is assessed by fitting the ELGM model
#' with and without the moderation (e.g. `V(t1) ~ z`) and comparing the marginal likelihoods with
#' [compare()] (a Bayes factor, which depends on the prior of the effect).
#' @param object A maximum likelihood fit of `birt(..., engine = "rtmb")` (`rtmb_control(method = "aghq")`).
#' @param p Quantile levels.
#' @param type `"conditional"` (covariates fixed) or `"unconditional"`.
#' @return A persons x effects matrix of influence functions (columns `target|predictor|p`) of
#'   class `birt_qscore`, with attribute `effects`: `est`, `se` (model-based), `se.robust`
#'   (sandwich).
#' @references
#' Debelak, R., Pawel, S., Strobl, C., & Merkle, E. C. (2022). Score-based measurement invariance
#' checks for Bayesian maximum-a-posteriori estimates in item response theory. *British Journal
#' of Mathematical and Statistical Psychology, 75*(3), 728-752. \doi{10.1111/bmsp.12275}
#'
#' Hampel, F. R. (1974). The influence curve and its role in robust estimation. *Journal of the
#' American Statistical Association, 69*, 383-393.
#'
#' Huber, P. J. (1967). The behavior of maximum likelihood estimates under nonstandard
#' conditions. *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and
#' Probability, 1*, 221-233.
#'
#' White, H. (1982). Maximum likelihood estimation of misspecified models. *Econometrica, 50*, 1-25.
#'
#' Cook, R. D. (1977). Detection of influential observation in linear regression.
#' *Technometrics, 19*, 15-18.
#'
#' Graubard, B. I., & Korn, E. L. (1999). Predictive margins with survey data. *Biometrics, 55*,
#' 652-659. \doi{10.1111/j.0006-341X.1999.00652.x}
#'
#' Zeileis, A., & Hornik, K. (2007). Generalized M-fluctuation tests for parameter instability.
#' *Statistica Neerlandica, 61*, 488-508. \doi{10.1111/j.1467-9574.2007.00371.x}
#' @export
quantile_score <- function(object, p = object$control$p, type = c("conditional", "unconditional")) {
  type <- match.arg(type)
  if (!inherits(object, "birt_rtmb")) stop("quantile_score() needs a fit of birt(engine = \"rtmb\")", call. = FALSE)
  need_ml(object, "quantile_score()")
  Sc <- estfun.birt_rtmb(object)                                    # N x parameters
  qf <- qeffect_fun(object$S, p)
  if (!qf$n) stop("the model has no quantile effects", call. = FALSE)
  d <- delta(object, qf$f, cov = TRUE); N <- nrow(Sc)
  keep <- !is.na(d$se)                                              # effects fixed by the model (e.g. a loading of -1) have no influence
  C <- Sc %*% object$rtmb$vcov %*% t(d$J[keep, , drop = FALSE])     # IF_j / N
  se <- d$se[keep]
  if (type == "unconditional") {
    H <- t(vapply(seq_len(N), function(j) unlist(qeffect_fun(object$S, p, only = j)$f(object$rtmb$par)), numeric(length(d$est))))[, keep, drop = FALSE]
    D <- sweep(H, 2, d$est[keep])                                   # h_j - beta: covariate sampling
    C <- C + D / N
    se <- sqrt(se^2 + colSums(D^2) / N^2)
  }
  IF <- N * C; colnames(IF) <- d$names[keep]; rownames(IF) <- NULL
  parts <- do.call(rbind, strsplit(colnames(IF), "|", fixed = TRUE))
  eff <- data.frame(target = parts[, 1], predictor = parts[, 2], p = as.numeric(parts[, 3]), est = d$est[keep],
                    se = se, se.robust = sqrt(colSums(C^2)), stringsAsFactors = FALSE)
  structure(IF, effects = eff, type = type, class = c("birt_qscore", "matrix"))
}

#' @rdname quantile_score
#' @param x A `birt_qscore` object.
#' @param ... Unused.
#' @export
print.birt_qscore <- function(x, ...) {
  cat(sprintf("Influence functions of quantile effects (%s): %d persons x %d effects\n", attr(x, "type"), nrow(x), ncol(x)))
  print_table(attr(x, "effects"), 4)
  invisible(x)
}

#' @rdname quantile_score
#' @param moderator One value per person: numeric, or a factor (logical and character are converted).
#' @param test `"auto"` (`"DM"` for a numeric moderator, `"LM"` for a factor), `"DM"`, `"CvM"`,
#'   `"maxLM"` or `"LM"`, as in [score_test()].
#' @details `score_test(qs, moderator)` orders each effect's influence function by the moderator
#'   and compares its cumulative sum, standardized by its variance, with a Brownian bridge. Under
#'   parameter invariance the partial sums of any fixed linear combination of the scores converge to
#'   a Brownian bridge (Zeileis & Hornik, 2007), so the test is valid; it detects changes of the
#'   parameters in the direction that moves the effect (to first order), one effect at a time, with
#'   Holm-adjusted p-values. It needs `type = "conditional"`: the covariate term of the
#'   unconditional influence function varies with the covariates by construction.
#' @export
score_test.birt_qscore <- function(object, moderator, test = c("auto", "DM", "CvM", "maxLM", "LM"), ...) {
  test <- match.arg(test)
  if (!requireNamespace("strucchange", quietly = TRUE)) stop("score_test() needs the 'strucchange' package", call. = FALSE)
  if (attr(object, "type") != "conditional") stop("score tests use quantile_score(type = \"conditional\")", call. = FALSE)
  N <- nrow(object)
  if (is.logical(moderator) || is.character(moderator)) moderator <- factor(moderator)
  if (length(moderator) != N) stop(sprintf("moderator has %d values for %d persons", length(moderator), N), call. = FALSE)
  if (anyNA(moderator)) stop("moderator has missing values", call. = FALSE)
  if (test == "auto") test <- if (is.factor(moderator)) "LM" else "DM"
  if (test == "LM" && !is.factor(moderator)) stop("test = \"LM\" is for a factor moderator", call. = FALSE)
  res <- t(vapply(seq_len(ncol(object)), function(k) bb_test(unclass(object)[, k, drop = FALSE], moderator, test), numeric(2)))
  out <- data.frame(parameters = colnames(object), n = 1L, statistic = res[, 1], p.value = res[, 2], row.names = NULL)
  out$p.holm <- stats::p.adjust(out$p.value, "holm")
  structure(out, class = c("birt_score_test", "data.frame"), test = test, decorrelate = "none",
            moderator = if (is.factor(moderator)) sprintf("factor with %d levels", nlevels(moderator)) else "numeric", clusters = NULL)
}
