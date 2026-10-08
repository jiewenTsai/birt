# Casewise scores of RTMB maximum likelihood fits and score-based tests of parameter invariance ---
#
# The score of person i is the derivative of that person's marginal log-likelihood (latent variables
# integrated out by the adaptive Gauss-Hermite rule of the fit, with the nodes of the final objective
# held fixed) with respect to the free parameters, on the internal (unconstrained) scale of the fit:
# log loadings for loadings constrained positive, log SDs, atanh correlations, ordinal thresholds as
# t1 and log increments. RTMB tapes the vector of per-person log-likelihoods and returns its Jacobian.

# one row per element of the free parameter vector: internal name, readable label, indicator, block
rtmb_par_info <- function(S, pf) {
  lens <- lengths(pf); out <- list()
  add <- function(par, i, label, item, block) out[[length(out) + 1]] <<- data.frame(par = sprintf("%s[%d]", par, i), label = label,
                                                                                    item = item, block = block, stringsAsFactors = FALSE)
  ldrow <- function(kind, i) {
    r <- S$ld[S$ld$kind == kind & S$ld$idx == i, ]
    list(label = if (nrow(r) == 1) sprintf("%s =~ %s", r$lhs, r$rhs) else sprintf("%s =~ (%s)", r$lhs[1], paste(r$rhs, collapse = ", ")),
         item = if (nrow(r) == 1) r$rhs else NA_character_, block = if (all(r$type %in% c("b", "o"))) "measurement" else "rt")
  }
  fs <- S$fs
  for (n in names(pf)) for (i in seq_len(lens[[n]])) {
    switch(n,
      d = add(n, i, sprintf("%s ~1", S$B[i]), S$B[i], "measurement"),
      lpos = , lreal = { r <- ldrow(if (n == "lpos") "pos" else "real", i); add(n, i, r$label, r$item, r$block) },
      xi = add(n, i, sprintf("%s ~1", S$C[i]), S$C[i], "rt"),
      lsig = add(n, i, sprintf("log sd(%s)", S$C[i]), S$C[i], "rt"),
      eps = add(n, i, sprintf("skew(%s)", S$C[S$sh[i]]), S$C[S$sh[i]], "rt"),
      ldl = add(n, i, sprintf("log tail(%s)", S$C[S$sh[i]]), S$C[S$sh[i]], "rt"),
      kap = add(n, i, sprintf("log sd(%s) ~ %s", S$C[S$kap$j[i]], S$kap$m[i]), S$C[S$kap$j[i]], "rt"),
      eta = add(n, i, sprintf("skew(%s) ~ %s", S$C[S$eta$j[i]], S$eta$m[i]), S$C[S$eta$j[i]], "rt"),
      beta1 = add(n, i, sprintf("%s ~ %s", fs[1], S$f[[1]]$covs[i]), NA, "structural"),
      beta2 = add(n, i, sprintf("%s ~ %s", fs[2], S$f[[2]]$covs[i]), NA, "structural"),
      gam = add(n, i, sprintf("%s ~ %s", fs[2], fs[1]), NA, "structural"),
      atr = add(n, i, sprintf("atanh cor(%s, %s)", fs[1], fs[2]), NA, "structural"),
      lsd1 = add(n, i, sprintf("log sd(%s)", fs[1]), NA, "structural"),
      lsd2 = add(n, i, sprintf("log sd(%s)", fs[2]), NA, "structural"),
      alpha1 = add(n, i, sprintf("log sd(%s) ~ %s", fs[1], S$f[[1]]$scale[i]), NA, "structural"),
      alpha2 = add(n, i, sprintf("log sd(%s) ~ %s", fs[2], S$f[[2]]$scale[i]), NA, "structural"),
      eps2 = add(n, i, sprintf("skew(%s)", fs[2]), NA, "structural"),
      ldl2 = add(n, i, sprintf("log tail(%s)", fs[2]), NA, "structural"),
      eta2 = add(n, i, sprintf("skew(%s) ~ %s", fs[2], S$f[[2]]$skew[i]), NA, "structural"),
      t1 = { j <- which(vapply(S$tix, function(x) i %in% x$t1, TRUE)); add(n, i, sprintf("b1(%s)", S$O[j]), S$O[j], "measurement") },
      linc = { j <- which(vapply(S$tix, function(x) i %in% x$linc, TRUE)); k <- match(i, S$tix[[j]]$linc) + 1; add(n, i, sprintf("log(b%d - b%d)(%s)", k, k - 1, S$O[j]), S$O[j], "measurement") },
      pb = { j <- which(vapply(S$tix, function(x) i %in% x$pb, TRUE)); k <- match(i, S$tix[[j]]$pb); add(n, i, sprintf("b%d(%s)", k, S$O[j]), S$O[j], "measurement") },
      lr = { j <- which(vapply(S$tix, function(x) i %in% x$lr, TRUE)); k <- match(i, S$tix[[j]]$lr) + 1; add(n, i, sprintf("log(a%d / a1)(%s)", k, S$O[j]), S$O[j], "measurement") },
      af = , bf = , ac = , bc = { m <- S$mod[S$mod$par == n & S$mod$idx == i, ]
        lab <- sprintf("E(%s) ~ %s", if (nrow(m) == 1) m$target else sprintf("(%s)", paste(m$target, collapse = ", ")),
                       if (m$kind[1] == "alpha") paste0(m$mod[1], ":", m$factor[1]) else m$mod[1])
        add(n, i, lab, if (nrow(m) == 1) m$target else NA_character_, if (all(m$target %in% c(S$B, S$O))) "measurement" else "rt") },
      add(n, i, sprintf("%s[%d]", n, i), NA, "structural"))
  }
  do.call(rbind, out)
}


# The score tools (estfun(), score_test(), quantile_score(), dif_tree()) are maximum likelihood
# methods: casewise scores at the ML estimate, which sum to zero, cumulated into a Brownian bridge
# (Zeileis & Hornik, 2007; Merkle & Zeileis, 2013). At the ELGM posterior mode the scores sum to
# minus the gradient of the log prior, and the reported ELGM estimates are posterior means, not
# the mode; MAP versions need centred scores and their error rates depend on the covariance used
# (Debelak, Pawel, Strobl & Merkle, 2022). So these tools take the ML fit of the same model.
need_ml <- function(x, what) {
  refit <- "refit the same model with control = rtmb_control(method = \"aghq\") (maximum likelihood) and use that fit"
  bayes <- "; the Bayesian route is to fit the ELGM model with and without the effect (e.g. E(y1) ~ z) and compare() them (Bayes factor from the marginal likelihoods)"
  if (!is.null(x$rtmb$elgm) || identical(x$control$method, "elgm"))
    stop(what, " is a maximum likelihood tool (casewise scores at the ML estimate) and this is an ELGM fit: ", refit, bayes, call. = FALSE)
  if (!identical(x$rtmb$method, "aghq"))
    stop(what, " needs the AGHQ maximum likelihood fit (its scores use the quadrature nodes of the final objective): ", refit, call. = FALSE)
  invisible(TRUE)
}

unpack_par <- function(x, pf) {
  out <- list(); o <- 0
  for (n in names(pf)) { L <- length(pf[[n]]); out[[n]] <- x[o + seq_len(L)]; o <- o + L }
  out
}

#' Casewise scores of an RTMB fit
#'
#' `estfun()` returns, for a maximum likelihood fit of [birt()] with `engine = "rtmb"` (adaptive
#' Gauss-Hermite quadrature, also for `ordered = TRUE`), the matrix of casewise scores: one row per
#' person, one column per free parameter, the derivative of that person's marginal log-likelihood
#' (latent variables integrated out, with the quadrature nodes of the final objective) at the
#' estimates. The columns sum to the negative gradient of the objective, about zero. Parameters are
#' on the internal scale of the fit (log of loadings constrained positive, log SDs, atanh of
#' correlations; for ordinal items log a, the first threshold and log increments); the attribute
#' `parameters` describes each column (`label`, indicator `item`, `block`: `"measurement"` for
#' binary and ordinal items, `"rt"` for continuous indicators, `"structural"`). These are the
#' estimating functions used by [score_test()] and by `strucchange::gefp()`. An ELGM fit is
#' refused (see "ELGM fits" in [score_test()]).
#' @param x A fit of `birt(..., engine = "rtmb")` (maximum likelihood, `rtmb_control(method = "aghq")`).
#' @param ... Unused.
#' @return A persons x parameters matrix.
#' @exportS3Method sandwich::estfun
estfun.birt_rtmb <- function(x, ...) {
  need_ml(x, "estfun()")
  S <- x$S; pf <- x$rtmb$par
  ND <- node_data(S, x$rtmb$centers, x$rtmb$k); N <- S$N; G <- ND$G
  info <- rtmb_par_info(S, pf)
  f <- function(v) {
    p <- unpack_par(v, pf)
    Fl <- lat_values(p, S, ND$z1, ND$z2, ND$rows)
    ll <- ind_ll(p, S, Fl, ND$rows) + ND$lw
    lacc <- ll[seq_len(N)]
    for (h in seq_len(G)[-1]) lacc <- RTMB::logspace_add(lacc, ll[(h - 1) * N + seq_len(N)])
    lacc + ND$cst
  }
  x0 <- unlist(x$rtmb$par, use.names = FALSE)
  tape <- RTMB::MakeTape(f, x0)
  J <- tape$jacobian(x0)
  colnames(J) <- info$label; rownames(J) <- NULL
  attr(J, "parameters") <- info
  J
}

#' @rdname estfun.birt_rtmb
#' @return `bread()`: the inverse of the average observed information (internal scale), as in the
#'   'sandwich' package.
#' @exportS3Method sandwich::bread
bread.birt_rtmb <- function(x, ...) {
  V <- x$rtmb$vcov; n <- x$data$N
  if (is.null(V) || anyNA(V)) stop("no covariance matrix in the fit", call. = FALSE)
  B <- V * n; nm <- colnames(estfun.birt_rtmb(x)); dimnames(B) <- list(nm, nm); B
}

#' Score-based tests of parameter invariance
#'
#' Tests whether item or structural parameters change along a person-level moderator (DIF,
#' measurement invariance, item position effects) from a single maximum likelihood fit
#' (`birt(..., engine = "rtmb")`, also with `ordered = TRUE`): the casewise scores
#' ([estfun.birt_rtmb()]) are ordered by the moderator, cumulated after decorrelation and compared
#' with a Brownian bridge (Merkle & Zeileis, 2013).
#'
#' **Moderators derived from the response times.** The marginal distribution of the continuous
#' indicators (response times) does not involve the parameters of the binary or ordinal items, so
#' in a joint model their scores have mean zero given any function of the times (for example the
#' total log time): the measurement parameters can be tested along such a moderator, provided the
#' scores are decorrelated within the tested block (`decorrelate = "block"`, the default). The
#' parameters of the continuous indicators cannot be tested along a moderator built from them.
#'
#' **DIF and impact.** A moderator related to the latent variable (impact) shifts all intercepts
#' and is flagged as instability. To test DIF beyond impact, add the moderator to the latent
#' regression (`f ~ z` in the syntax, or `E(f) ~ z` for ordinal items) and test that fit;
#' `score_test()` gives a message when the moderator correlates with the factor scores and is not a
#' covariate of the model.
#'
#' **Clustered data.** With `cluster`, the scores are summed within clusters, which is exact for a
#' moderator that is constant within clusters; a person-level moderator in clustered data is refused.
#'
#'
#' **ELGM fits.** `score_test()` is a maximum likelihood tool and refuses an ELGM fit (the default of
#' [birt()]): refit the same model with `control = rtmb_control(method = "aghq")`. The tests
#' need casewise scores that sum to zero at the estimate; at the ELGM posterior mode they sum to
#' minus the gradient of the log prior, the reported ELGM estimates are posterior means rather
#' than the mode, and MAP versions of the score tests need centred scores, with error rates that
#' depend on the covariance estimate (Debelak, Pawel, Strobl & Merkle, 2022). The Bayesian
#' route for one effect is to fit the ELGM model with and without it (e.g. `E(y1) ~ z` for
#' uniform DIF of `y1` along a numeric `z`) and compare the marginal likelihoods with
#' [compare()] (a Bayes factor, which depends on the prior of the effect, `dpriors(beta = )`).
#'
#' birtRcpp has a function of the same name for its models; with both packages attached, the one
#' attached last is found first, and `birt::score_test()` dispatches to birtRcpp for its models.
#' @references
#' Merkle, E. C., & Zeileis, A. (2013). Tests of measurement invariance without subgroups: A
#' generalization of classical methods. *Psychometrika, 78*(1), 59-82. \doi{10.1007/s11336-012-9302-4}
#'
#' Debelak, R., Pawel, S., Strobl, C., & Merkle, E. C. (2022). Score-based measurement invariance
#' checks for Bayesian maximum-a-posteriori estimates in item response theory. *British Journal
#' of Mathematical and Statistical Psychology, 75*(3), 728-752. \doi{10.1111/bmsp.12275}
#' @param object A maximum likelihood fit of `birt(..., engine = "rtmb")` (`rtmb_control(method = "aghq")`).
#' @param moderator One value per person: numeric (continuous or ordered) or a factor (logical
#'   and character are turned into factors).
#' @param parm Parameters to test: `"measurement"` (loadings and intercepts or
#'   thresholds of the binary and ordinal items, plus their moderation effects; default), `"rt"`
#'   (parameters of the continuous indicators, including cross-loadings on them), `"structural"`, a
#'   regular expression on the parameter labels, or column indices.
#' @param by_item `TRUE`: also test the parameters of each indicator separately, with Holm-adjusted
#'   p-values.
#' @param test `"auto"` (`"DM"` for numeric moderators, `"LM"` for factors), `"DM"`
#'   (double maximum), `"CvM"` (Cramér-von Mises), `"maxLM"` (supremum LM, 10\% trimming) or `"LM"`.
#' @param decorrelate `"block"` (within the tested parameters) or `"full"` (all parameters, then the
#'   tested components; not valid for moderators built from the response times).
#' @param cluster Optional cluster identifiers (see Details).
#' @param impact_check Message when the moderator correlates with the factor scores but is not a
#'   covariate of the model.
#' @param ... Passed to the method.
#' @return A data frame of class `birt_score_test`: `parameters`, `n`, `statistic`, `p.value` and,
#'   with `by_item`, `p.holm`.
#' @examples
#' \donttest{
#' d <- sim_rtirt(N = 400, K = 6)
#' f <- birt(rtirt_syntax(paste0("y", 1:6), paste0("t", 1:6)), d,
#'           control = rtmb_control(method = "aghq"), progress = FALSE)
#' score_test(f, d$x1 > 0, by_item = TRUE)            # a moderator unrelated to the items
#' score_test(f, rowSums(d[paste0("t", 1:6)]))        # along the total log time
#' }
#' @export
score_test <- function(object, ...) UseMethod("score_test")

#' @export
score_test.default <- function(object, ...) {
  if (inherits(object, "rtirt") && isNamespaceLoaded("birtRcpp"))          # a birtRcpp model: its own function
    return(getExportedValue("birtRcpp", "score_test")(object, ...))
  stop("score_test() needs a maximum likelihood fit of birt(..., engine = \"rtmb\")", call. = FALSE)
}

#' @rdname score_test
#' @export
score_test.birt_rtmb <- function(object, moderator, parm = "measurement", by_item = FALSE,
                                 test = c("auto", "DM", "CvM", "maxLM", "LM"), decorrelate = c("block", "full"),
                                 cluster = NULL, impact_check = TRUE, ...) {
  test <- match.arg(test); decorrelate <- match.arg(decorrelate)
  if (!requireNamespace("strucchange", quietly = TRUE)) stop("score_test() needs the 'strucchange' package", call. = FALSE)
  need_ml(object, "score_test()")
  S <- estfun.birt_rtmb(object); info <- attr(S, "parameters"); N <- nrow(S); nm <- colnames(S)
  if (is.logical(moderator) || is.character(moderator)) moderator <- factor(moderator)
  if (length(moderator) != N) stop(sprintf("moderator has %d values for %d persons", length(moderator), N), call. = FALSE)
  if (anyNA(moderator)) stop("moderator has missing values; drop those persons before fitting", call. = FALSE)
  if (impact_check) {
    zz <- as.numeric(if (is.factor(moderator)) moderator != levels(moderator)[1] else moderator)
    X <- object$S$X
    explained <- length(X) && ncol(X) && max(abs(stats::cor(zz, X))) > 0.99
    fs <- setdiff(names(object$scores), c("id", grep("_psd$", names(object$scores), value = TRUE)))
    for (fn in fs) {
      r <- stats::cor(object$scores[[fn]], zz)
      if (!explained && abs(r) > 0.1 && stats::cor.test(object$scores[[fn]], zz)$p.value < 0.001)
        message(sprintf(paste0("the moderator correlates %.2f with the scores of %s: differences in its mean (impact) shift every ",
                               "intercept; to test DIF beyond impact, add the moderator to the latent regression (%s ~ z, or E(%s) ~ z) ",
                               "and test that fit"), r, fn, fn, fn))
    }
  }
  if (!is.null(cluster)) {
    if (length(cluster) != N) stop("cluster needs one value per person", call. = FALSE)
    same <- tapply(seq_len(N), cluster, function(i) length(unique(as.character(moderator[i]))) == 1)
    if (!all(same)) stop("with cluster, the moderator must be constant within clusters (a person-level moderator in clustered data is not covered)", call. = FALSE)
    first <- !duplicated(cluster)
    S <- rowsum(S, cluster, reorder = FALSE); moderator <- moderator[first]; N <- nrow(S)
  }
  idx <- if (is.numeric(parm)) parm else switch(parm,
    measurement = which(info$block == "measurement"), rt = which(info$block == "rt"),
    structural = which(info$block == "structural"), grep(parm, nm))
  if (!length(idx)) stop("no parameters match parm", call. = FALSE)
  if (test == "auto") test <- if (is.factor(moderator)) "LM" else "DM"
  if (test == "LM" && !is.factor(moderator)) stop("test = \"LM\" is for a factor moderator", call. = FALSE)
  one <- function(j) if (decorrelate == "block") bb_test(S[, j, drop = FALSE], moderator, test) else bb_test(S, moderator, test, parm = j)
  out <- data.frame(parameters = if (is.character(parm) && length(parm) == 1) parm else "selected", n = length(idx), t(one(idx)), row.names = NULL)
  if (by_item) {
    items <- unique(stats::na.omit(info$item[idx]))
    per <- lapply(items, function(it) { j <- idx[which(info$item[idx] == it)]; c(it, length(j), one(j)) })
    if (length(per)) {
      per <- do.call(rbind, per)
      pt <- data.frame(parameters = per[, 1], n = as.integer(per[, 2]), statistic = as.numeric(per[, 3]), p.value = as.numeric(per[, 4]))
      pt$p.holm <- stats::p.adjust(pt$p.value, "holm"); out$p.holm <- NA
      out <- rbind(out, pt)
    }
  }
  structure(out, class = c("birt_score_test", "data.frame"), test = test, decorrelate = decorrelate,
            moderator = if (is.factor(moderator)) sprintf("factor with %d levels", nlevels(moderator)) else "numeric",
            clusters = if (!is.null(cluster)) N)
}
# cumulated scores ordered by the moderator against a Brownian bridge (Zeileis & Hornik, 2007):
# statistic and p-value for the columns parm of S (decorrelated with all columns of S)
bb_test <- function(S, moderator, test, parm = seq_len(ncol(S))) {
  obj <- structure(list(S = S), class = "birt_scores")
  gp <- strucchange::gefp(obj, fit = NULL, scores = function(x, ...) x$S, order.by = moderator, parm = parm, sandwich = FALSE)
  fun <- switch(test, DM = strucchange::maxBB, CvM = strucchange::meanL2BB, maxLM = strucchange::supLM(0.1),
                LM = strucchange::catL2BB(gp))
  r <- suppressWarnings(strucchange::sctest(gp, functional = fun))
  c(statistic = unname(r$statistic), p.value = unname(r$p.value))
}

#' @export
coef.birt_scores <- function(object, ...) stats::setNames(rep(0, ncol(object$S)), colnames(object$S))

#' @export
print.birt_score_test <- function(x, digits = 3, ...) {
  cat(sprintf("Score-based test of parameter invariance (%s test, %s moderator, %s decorrelation%s)\n", attr(x, "test"),
              attr(x, "moderator"), attr(x, "decorrelate"), if (!is.null(attr(x, "clusters"))) sprintf(", %d clusters", attr(x, "clusters")) else ""))
  d <- as.data.frame(unclass(x)); attributes(d)[c("test", "decorrelate", "moderator", "clusters")] <- NULL
  d$statistic <- round(d$statistic, digits); d$p.value <- format.pval(d$p.value, digits = digits, eps = 1e-4)
  if (!is.null(d$p.holm)) d$p.holm <- ifelse(is.na(d$p.holm), "", format.pval(d$p.holm, digits = digits, eps = 1e-4))
  print(d, row.names = FALSE)
  invisible(x)
}
