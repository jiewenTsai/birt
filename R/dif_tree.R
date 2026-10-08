# DIF trees: model-based recursive partitioning of a birt model --------------------------------------

#' DIF tree
#'
#' Recursively splits the persons by the partitioning variables where the measurement parameters of
#' the binary or ordinal items (loadings, intercepts, thresholds) are not invariant: model-based
#' recursive partitioning (Zeileis, Hothorn & Hornik, 2008; for IRT, Strobl, Kopf & Zeileis, 2015).
#' At each node the model is fitted by maximum likelihood (`engine = "rtmb"`) and the casewise scores
#' of the measurement parameters ([estfun.birt_rtmb()]) are tested along every partitioning variable
#' with fluctuation tests, as in [score_test()]; the most significant variable (Bonferroni-adjusted) is
#' split if it is below `alpha`.
#'
#' **Impact.** With `impact = TRUE` (default) the partitioning variables are covariates of every latent
#' variable in each node (`f ~ z`, centred and scaled within the node; factors as dummies), so
#' differences in the latent means between groups are modelled and the tree looks for DIF beyond
#' impact. With `impact = FALSE` any variable related to the latent variables splits the tree.
#'
#' **Tests.** Numeric partitioning variables are tested by the supLM (maxLM) statistic, ordered
#' factors (including the quantile groups of `nbins`) by its version for ordinal variables
#' (maxLM_o, `mob_control(ordinal = "L2")`; Merkle, Fan & Zeileis, 2014) and unordered factors by
#' the chi-square statistic.
#'
#' **ELGM.** The tree is a maximum likelihood method: each node is fitted with
#' `rtmb_control(method = "aghq")` (the default `control`) and an ELGM `control` is refused, since
#' the fluctuation tests need casewise scores at the ML estimate (see "ELGM fits" in
#' [score_test()]). There is no ELGM counterpart of the recursive search; a single candidate split
#' can be checked by comparing ELGM fits with and without `E(y) ~ z` ([compare()]).
#'
#' Only the scores of the measurement block are tested (decorrelated within the block), so
#' partitioning variables built from the response times (e.g. the total log time) are allowed.
#' Parameters are on the scale of each node, so compare nodes through the tests rather than raw
#' parameter differences.
#' @param model Model syntax, as in [birt()].
#' @param data A data frame with the indicators.
#' @param partition A data frame of partitioning variables (one row per person; numeric, ordered or
#'   factor).
#' @param impact Model latent mean differences along the partitioning variables.
#' @param nbins Numeric partitioning variables with more than `nbins` distinct values are cut into
#'   `nbins` groups at their quantiles (an ordered factor; the tree prints the groups as intervals).
#'   A split is then searched among the `nbins - 1` group boundaries (two maximum likelihood fits
#'   each) instead of every observed value (about two fits per person in the node), which makes
#'   the tree minutes faster; the cost is that a cut point is located only up to its quantile group.
#'   `NULL` keeps them numeric (the supLM test and the exhaustive search of [partykit::mob()]).
#' @param alpha,minsize,maxdepth Significance level (Bonferroni over variables), minimum node size
#'   and maximum depth (see [partykit::mob_control()]).
#' @param ... Passed to [birt()] (e.g. `ordered`, `family`; `control` defaults to the maximum likelihood fit,
#'   `rtmb_control(method = "aghq")`, whose scores the tests use).
#' @return A `modelparty` object (partykit): `print()`, `plot()`, `coef()` (measurement parameters per
#'   terminal node, internal scale), `strucchange::sctest(tree, node = )` for the tests at a node.
#' @references
#' Debelak, R., Pawel, S., Strobl, C., & Merkle, E. C. (2022). Score-based measurement invariance
#' checks for Bayesian maximum-a-posteriori estimates in item response theory. *British Journal
#' of Mathematical and Statistical Psychology, 75*(3), 728-752. \doi{10.1111/bmsp.12275}
#'
#' Zeileis, A., Hothorn, T., & Hornik, K. (2008). Model-based recursive partitioning. *Journal of
#' Computational and Graphical Statistics, 17*, 492-514. \doi{10.1198/106186008X319331}
#'
#' Merkle, E. C., Fan, J., & Zeileis, A. (2014). Testing for measurement invariance with respect to
#' an ordinal variable. *Psychometrika, 79*, 569-584. \doi{10.1007/s11336-013-9376-7}
#'
#' Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method for detecting differential
#' item functioning in the Rasch model. *Psychometrika, 80*, 289-316. \doi{10.1007/s11336-013-9388-3}
#' @examples
#' \donttest{
#' if (requireNamespace("partykit", quietly = TRUE) && requireNamespace("RTMB", quietly = TRUE)) {
#'   d <- sim_rtirt(N = 600, K = 6, seed = 1)
#'   z <- data.frame(group = factor(sample(c("A", "B"), 600, TRUE)), x = rnorm(600))
#'   dif_tree(rtirt_syntax(paste0("y", 1:6), paste0("t", 1:6)), d, z)   # no DIF: a single node
#' }
#' }
#' @export
dif_tree <- function(model, data, partition, impact = TRUE, nbins = 10, alpha = 0.05, minsize = NULL, maxdepth = Inf, ...) {
  for (pk in c("partykit", "RTMB")) if (!requireNamespace(pk, quietly = TRUE)) stop("dif_tree() needs the '", pk, "' package", call. = FALSE)
  partition <- as.data.frame(partition); data <- as.data.frame(data)
  if (nrow(partition) != nrow(data)) stop(sprintf("partition has %d rows for %d persons", nrow(partition), nrow(data)), call. = FALSE)
  if (anyNA(partition)) stop("partition has missing values", call. = FALSE)
  for (v in names(partition)) if (is.character(partition[[v]]) || is.logical(partition[[v]])) partition[[v]] <- factor(partition[[v]])
  if (!is.null(nbins) && (!is.numeric(nbins) || length(nbins) != 1 || nbins < 2)) stop("nbins must be NULL or a number >= 2", call. = FALSE)
  if (!is.null(nbins)) for (v in names(partition))
    if (is.numeric(partition[[v]]) && length(unique(partition[[v]])) > nbins) partition[[v]] <- quantile_groups(partition[[v]], nbins)
  bargs <- list(...)                                                       # arguments of birt()
  bargs$control <- bargs$control %||% rtmb_control(method = "aghq")        # the scores need the maximum likelihood fit
  if (!identical(bargs$control$method, "aghq"))
    stop("dif_tree() fits each node by maximum likelihood and tests its casewise scores: give control = rtmb_control(method = \"aghq\", ...) (the default) or leave control out; ",
         "an ELGM tree has no counterpart in birt (for one split, compare() ELGM fits with and without E(y) ~ z)", call. = FALSE)
  sp <- build_spec(model, data, bargs$binary, bargs$family, bargs$ordered, bargs$itemtype %||% "grm")
  if (!length(c(sp$binary, sp$ordinal))) stop("dif_tree() tests binary or ordinal items; the model has none", call. = FALSE)
  vars <- unique(c(sp$binary, sp$ordinal, sp$cont, sp$covs, sp$modvars))
  clash <- intersect(names(partition), vars)
  if (length(clash)) stop("partitioning variables are also variables of the model: ", paste(clash, collapse = ", "), call. = FALSE)
  fit_node <- function(y, x, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
    dn <- as.data.frame(matrix(as.matrix(y), ncol = length(vars), dimnames = list(NULL, vars)))
    X <- as.matrix(x)[, -1, drop = FALSE]
    X <- X[, apply(X, 2, stats::sd) > 0, drop = FALSE]                       # constant within the node after a split
    m <- model
    if (impact && ncol(X)) {
      X <- scale(X); colnames(X) <- paste0(".z", seq_len(ncol(X))); dn <- cbind(dn, X)
      m <- paste(c(model, sprintf("%s ~ %s", sp$factors, paste(colnames(X), collapse = " + "))), collapse = "\n")
    }
    f <- suppressWarnings(suppressMessages(do.call(birt, c(list(m, dn, engine = "rtmb", progress = FALSE), bargs))))
    S <- estfun.birt_rtmb(f); j <- which(attr(S, "parameters")$block == "measurement")
    cf <- unlist(f$rtmb$par, use.names = FALSE)[j]; names(cf) <- colnames(S)[j]
    list(coefficients = cf, objfun = f$rtmb$objective, estfun = if (estfun) S[, j, drop = FALSE], object = if (object) f)
  }
  d <- partition; d$.y <- as.matrix(data[vars])
  pv <- paste(sprintf("`%s`", names(partition)), collapse = " + ")
  fo <- stats::as.formula(sprintf(".y ~ %s | %s", if (impact) pv else "1", pv))
  K <- length(c(sp$binary, sp$ordinal))
  ctrl <- partykit::mob_control(alpha = alpha, bonferroni = TRUE, minsize = minsize %||% max(100, 20 * K), maxdepth = maxdepth,
                                ordinal = "L2")
  partykit::mob(fo, data = d, fit = fit_node, control = ctrl)
}

# an ordered factor of nbins groups at the quantiles of x (fewer when quantiles tie)
quantile_groups <- function(x, nbins) {
  br <- unique(stats::quantile(x, seq(0, 1, length.out = nbins + 1), names = FALSE))
  cut(x, br, include.lowest = TRUE, ordered_result = TRUE, dig.lab = 3)
}
