#' The `birt` object and its accessors
#'
#' A JAGS fit of [birt()] is a list of class `c("birt", "birt_fit")` with these elements
#' and the accessors `scores()`, `reliability()`, `fit_indices()`, `pointwise_loglik()`,
#' `posterior_draws()`, `convergence()`, `jags_code()`, `compare()` (RTMB fits, including
#' [hgrm()], are of class `birt_rtmb`):
#' \describe{
#' \item{`estimates`}{lavaan-like parameter table: `lhs`, `op`, `rhs`, `family`,
#'   `prior`, posterior `est` (mean), `sd`, `q025`, `q975`, `rhat`, `ess`, ALD `scale`,
#'   spike-and-slab `p_incl` and `BF10` (posterior over prior odds of inclusion; the prior
#'   odds from `ssp_p`, or from the posterior mean of the learned inclusion probability of
#'   the parameter's group). Rows with `op` `~~` report variances (for ALD:
#'   the implied variance, with the scale in `scale`).}
#' \item{`scores`}{one row per person: EAP (posterior mean) and posterior SD of each
#'   latent variable.}
#' \item{`reliability`}{empirical reliability of the EAP scores, var(EAP) /
#'   (var(EAP) + mean PSD^2).}
#' \item{`items`}{item tables: binary items (a, b = -d/a, d) and continuous indicators
#'   (intercept, loadings, residual scale and variance, inclusion probability).}
#' \item{`loglik`}{draws x persons matrix of the marginal log-likelihood (latent
#'   variables integrated out); `loo`, `waic` the corresponding 'loo' objects;
#'   `fit_indices` Deviance, pD, DIC, WAIC, LOOIC.}
#' \item{`draws`}{the posterior draws ([coda::mcmc.list]), including person parameters.}
#' \item{`code`, `file`, `partab`, `spec`, `data`, `mcmc`, `secs`}{JAGS code and file,
#'   internal parameter table, parsed model, JAGS data, MCMC settings, run time.}
#' }
#' @name birt-class
#' @param object,x A `birt` (or `birt_quantile`) object.
#' @param ... Unused.
NULL

#' @rdname birt-class
#' @export
estimates <- function(object, ...) UseMethod("estimates")
#' @rdname birt-class
#' @export
estimates.birt <- function(object, ...) {
  e <- object$estimates
  cols <- c("lhs", "op", "rhs", "family", "prior", "est", "sd", "q025", "q975", "rhat", "ess", "scale", "p_incl", "BF10")
  if (any(!is.na(e$sd_adj))) cols <- c(cols, "sd_adj", "q025_adj", "q975_adj")   # ALD: Yang et al. (2016)
  e[, intersect(cols, names(e))]
}

#' @rdname birt-class
#' @export
scores <- function(object, ...) UseMethod("scores")
#' @rdname birt-class
#' @export
scores.birt_fit <- function(object, ...) object$scores

#' @rdname birt-class
#' @export
reliability <- function(object, ...) UseMethod("reliability")
#' @rdname birt-class
#' @export
reliability.birt_fit <- function(object, ...) object$reliability


#' @rdname birt-class
#' @export
fit_indices <- function(object, ...) UseMethod("fit_indices")
#' @rdname birt-class
#' @export
fit_indices.birt_fit <- function(object, ...) {
  if (is.null(object$fit_indices)) stop("no fit indices stored; refit with fit_indices = TRUE or use add_fit_indices()")
  object$fit_indices
}

#' @rdname birt-class
#' @export
pointwise_loglik <- function(object, ...) object$loglik

#' @rdname birt-class
#' @param persons Include the person parameters.
#' @export
posterior_draws <- function(object, persons = FALSE, ...) {
  M <- as.matrix(object$draws)
  if (!persons) M <- M[, !grepl(sprintf("^(%s)\\[", paste(latent_names(object), collapse = "|")), colnames(M)), drop = FALSE]
  M
}

#' @exportS3Method coda::as.mcmc
as.mcmc.birt_fit <- function(x, ...) x$draws

latent_names <- function(object) object$spec$factors %||% object$spec$factor

#' JAGS code of a model
#'
#' @param x A fitted `birt` object, or model syntax for [birt()].
#' @param data,binary,family,dp,ssp_p,ordered,itemtype As in [birt()] (when `x` is syntax).
#' @return The JAGS model code (character), printed with `cat()`.
#' @export
jags_code <- function(x, data = NULL, binary = NULL, family = NULL, dp = dpriors(), ssp_p = NULL, ordered = NULL, itemtype = "grm") {
  code <- if (inherits(x, "birt_fit")) x$code else {
    if (is.null(data)) stop("jags_code(): `data` is needed to build the code from syntax")
    sp <- build_spec(x, data, binary, family, ordered, itemtype); sp$dp <- dp; sp$ssp_p <- ssp_p
    build_code(sp)$code
  }
  cat(code)
  invisible(code)
}

#' @rdname birt-class
#' @export
coef.birt <- function(object, ...) {
  e <- object$estimates[!object$estimates$fixed, ]
  stats::setNames(e$est, paste0(e$lhs, e$op, e$rhs))
}

# the lines of the printed header that describe the model (both engines)
model_lines <- function(x) {
  sp <- x$spec
  fams <- c(vapply(sp$factors, function(f) sprintf("%s: %s", f, fam_label(sp$family$factor[[f]])), ""),
            if (length(sp$cont)) {
              tab <- table(vapply(sp$family$ind, fam_label, ""))
              sprintf("%d continuous: %s", length(sp$cont), paste(sprintf("%s (%d)", names(tab), tab), collapse = ", "))
            },
            if (length(sp$binary)) sprintf("%d binary: logit", length(sp$binary)),
            if (length(sp$ordinal)) ord_desc(sp))
  c(sprintf("  %d persons; latent variables: %s", x$data$N, paste(sp$factors, collapse = ", ")),
    sprintf("  Distributions: %s", paste(fams, collapse = "; ")),
    if (length(sp$covs)) sprintf("  Covariates: %s", paste(sp$covs, collapse = ", ")),
    if (length(sp$modvars)) sprintf("  Moderators: %s (%s)", paste(sp$modvars, collapse = ", "), moderation_txt(sp)))
}

header_lines <- function(x) {
  m <- x$mcmc
  nd <- coda::niter(x$draws) * coda::nchain(x$draws)
  c(sprintf("birt %s%s", pkg_version(), if (!is.null(x$label)) sprintf(" [%s]", x$label) else ""),
    model_lines(x),
    sprintf("  JAGS, conditional likelihood (persons sampled with the parameters): %d chains x %d iterations, burn-in %d, thin %d; %d draws; %.0f s",
            m$n_chains, m$n_iter, m$n_burn, m$thin, nd, x$secs))
}

#' @rdname birt-class
#' @export
print.birt <- function(x, ...) {
  cat(header_lines(x), sep = "\n")
  cat(sprintf("  Reliability: %s\n", paste(sprintf("%s %.3f", names(x$reliability), x$reliability), collapse = ", ")))
  if (!is.null(x$fit_indices))
    cat(sprintf("  LOOIC %.1f (SE %.1f), WAIC %.1f, DIC %.1f\n", x$fit_indices$LOOIC, x$fit_indices$SE_LOOIC,
                x$fit_indices$WAIC, x$fit_indices$DIC))
  mr <- max_na(x$summary_nodes$rhat)
  if (is.finite(mr) && mr > 1.05) cat(sprintf("  Warning: max Rhat %.2f > 1.05; see convergence()\n", mr))
  cat("Use summary() for the parameter table.\n")
  invisible(x)
}

# evidence category of a Bayes factor BF10 (Jeffreys' bands 1/10, 1/3, 3, 10, as in Lee &
# Wagenmakers, 2013); BF10 = Inf (p_incl = 1 in every draw) is strong evidence for nonzero
bf_evidence <- function(bf)
  as.character(cut(bf, c(-Inf, 1 / 10, 1 / 3, 3, 10, Inf), right = FALSE, include.lowest = TRUE,
                   labels = c("strong for 0", "moderate for 0", "inconclusive", "moderate for nonzero", "strong for nonzero")))

#' @rdname birt-class
#' @param digits Number of decimals.
#' @param irt Also show binary items in the IRT parameterization (a, b).
#' @export
summary.birt <- function(object, digits = 3, irt = TRUE, ...) {
  fit <- object; e <- fit$estimates; sp <- fit$spec; M <- as.matrix(fit$draws)
  bayes_cols <- e
  par <- partable(bayes_cols, TRUE, "posterior mean, SD and 95% interval; prior of each parameter", spec = sp)
  # IRT parameterization of the binary items (a = loading, b = -d / a), from the draws
  irt_rows <- NULL
  if (irt && length(sp$binary)) {
    f1 <- sp$factors[1]
    for (v in sp$binary) {
      ia <- which(e$op == "=~" & e$lhs == f1 & e$rhs == v); id <- which(e$op == "~1" & e$lhs == v)
      if (!length(ia) || !length(id) || is.na(e$node[id])) next
      a <- if (is.na(e$node[ia])) rep(e$est[ia], nrow(M)) else M[, e$node[ia]]; b <- -M[, e$node[id]] / a
      irt_rows <- rbind(irt_rows, data.frame(lhs = v, op = c("irt_a", "irt_b"), rhs = "", est = c(mean(a), mean(b)),
                                             se = c(stats::sd(a), stats::sd(b)), stringsAsFactors = FALSE))
    }
  }
  w <- item_wide(rbind(as.data.frame(unclass(par))[, c("lhs", "op", "rhs", "est", "se")], irt_rows), sp$ovs, sp$factors, sp$binary, sp$ordinal)
  fi <- fit$fit_indices
  fit_tab <- data.frame(N = fit$data$N, max_rhat = max_na(fit$summary_nodes$rhat))
  if (is.data.frame(fi) && nrow(fi)) fit_tab <- cbind(fit_tab, fi[, intersect(c("DIC", "pD", "WAIC", "LOOIC", "SE_LOOIC", "k_gt_0.7"), names(fi)), drop = FALSE])
  sel <- NULL; ss <- e[!is.na(e$p_incl), ]
  if (nrow(ss)) {
    if (is.null(ss$group)) ss$group <- fit$partab$group[match(ss$node, fit$partab$node)]
    nd <- coda::niter(fit$draws) * coda::nchain(fit$draws)
    pri <- vapply(seq_len(nrow(ss)), function(i) { o <- prior_odds(ss[i, ], M, sp); o / (1 + o) }, 0)
    ev <- bf_evidence(ss$BF10)
    sel <- btable(data.frame(parameter = paste(ss$lhs, ss$op, ss$rhs), est = ss$est, prior_incl = pri, p_incl = ss$p_incl,
                             BF10 = ss$BF10, selected = ifelse(ss$p_incl > 0.5, "nonzero", "zero"), evidence = ev),
                  sprintf("spike-and-slab: p_incl = P(not 0); selected by p_incl > .5; BF10 = posterior / prior odds%s",
                          if (is.null(sp$ssp_p)) "; prior_incl = posterior mean of the learned inclusion probability of the group (beta(1, 1) prior)" else sprintf(" (prior inclusion %s)", format(sp$ssp_p))))
  }
  structure(list(header = header_lines(fit),
                 fit = btable(fit_tab, "fit indices from the marginal likelihood (persons integrated out)"),
                 items = btable(w$est, sprintf("item parameters (posterior means)%s", if (irt && length(sp$binary)) "; binary items also as IRT a, b = -d / a" else "")),
                 se = btable(w$se, "posterior SDs of the item parameters"),
                 moderation = mod_summary(fit), precision = btable(fit$precision, "precision of the scores by tercile of the first moderator"),
                 parameters = par, selection = sel,
                 reliability = btable(as.data.frame(as.list(fit$reliability)), "reliability of the EAP scores: var(EAP) / (var(EAP) + mean PSD^2)")),
            class = c("summary.birt", "birt_summary"))
}

#' Convergence summary
#'
#' JAGS fits: R-hat and effective sample sizes, and the parameters with R-hat above
#' `threshold`. RTMB fits: optimizer code, maximum gradient and the AGHQ rounds.
#' @param object A `birt` or `birt_rtmb` fit.
#' @param threshold Rhat threshold for listing parameters.
#' @param ... Unused.
#' @return A data frame of parameters with Rhat above `threshold` (invisibly).
#' @export
convergence <- function(object, ...) UseMethod("convergence")

#' @rdname convergence
#' @export
convergence.birt_fit <- function(object, threshold = 1.05, ...) {
  sn <- object$summary_nodes
  print_table(data.frame(parameters = nrow(sn), max_rhat = max_na(sn$rhat),
                       n_rhat_gt_1.01 = sum(sn$rhat > 1.01, na.rm = TRUE),
                       n_rhat_gt_threshold = sum(sn$rhat > threshold, na.rm = TRUE),
                       min_ess = min(sn$ess[sn$ess > 0], na.rm = TRUE)))
  bad <- sn[!is.na(sn$rhat) & sn$rhat > threshold, ]
  if (nrow(bad)) print_table(bad[, c("node", "mean", "rhat", "ess")])
  invisible(bad)
}

#' Compare fitted models
#'
#' @param ... `birt` objects (named arguments give the model names), or one named list.
#'   For RTMB fits (`birt_rtmb`): -2 log L, number of parameters, AIC, BIC and their
#'   differences to the best model.
#' @param digits Number of decimals.
#' @return A data frame with reliability, DIC, WAIC, LOOIC and the ELPD difference to the
#'   best model (from [loo::loo_compare()]), invisibly.
#' @export
compare <- function(..., digits = 3) {
  fits <- list(...)
  if (length(fits) == 1 && is.list(fits[[1]]) && !inherits(fits[[1]], c("birt_fit", "birt_rtmb"))) fits <- fits[[1]]
  if (all(vapply(fits, inherits, logical(1), "birt_rtmb"))) return(compare_rtmb(fits, digits))
  if (any(vapply(fits, inherits, logical(1), "birt_rtmb"))) stop("compare() takes either JAGS fits or RTMB fits, not both")
  if (!all(vapply(fits, inherits, logical(1), "birt_fit"))) stop("compare() takes fitted birt or hgrm objects")
  if (is.null(names(fits)) || any(names(fits) == "")) names(fits) <- sprintf("model%d", seq_along(fits))
  fs <- unique(unlist(lapply(fits, function(f) names(f$reliability))))
  tab <- do.call(rbind, lapply(names(fits), function(m) {
    f <- fits[[m]]; fi <- f$fit_indices
    rel <- stats::setNames(as.list(f$reliability[fs]), paste0("rel_", fs))
    data.frame(model = m, rel, DIC = fi$DIC %||% NA, WAIC = fi$WAIC %||% NA, LOOIC = fi$LOOIC %||% NA,
               SE_LOOIC = fi$SE_LOOIC %||% NA, max_rhat = max_na(f$summary_nodes$rhat), check.names = FALSE)
  }))
  ok <- vapply(fits, function(f) !is.null(f$loo), logical(1))
  same_n <- length(unique(vapply(fits, function(f) f$data$N, 1))) == 1
  same_ind <- length(unique(vapply(fits, function(f) paste(sort(f$spec$ovs), collapse = ","), ""))) == 1
  if (sum(ok) > 1 && same_n && same_ind) {
    lc <- loo::loo_compare(lapply(fits[ok], `[[`, "loo"))
    tab$elpd_diff <- lc[match(tab$model, rownames(lc)), "elpd_diff"]
    tab$se_diff <- lc[match(tab$model, rownames(lc)), "se_diff"]
  }
  print_table(tab, digits)
  if (!same_ind) cat("Note: the models have different indicators; their criteria are not comparable.\n")
  invisible(tab)
}

compare_rtmb <- function(fits, digits) {
  if (is.null(names(fits)) || any(names(fits) == "")) names(fits) <- sprintf("model%d", seq_along(fits))
  bayes <- vapply(fits, function(f) !is.null(f$rtmb$elgm), TRUE)
  if (any(bayes) && !all(bayes)) stop("compare(): fits with method = \"elgm\" are compared by marginal likelihoods; refit all with the same method")
  if (all(bayes)) {
    tab <- do.call(rbind, lapply(names(fits), function(m) cbind(model = m, fits[[m]]$fit_indices[, c("method", "npar", "logML")])))
    tab$dlogML <- tab$logML - max(tab$logML)
    tab$BF_best <- formatC(exp(-tab$dlogML), digits = 3, format = "g")      # e.g. 6.26e+11
    if (length(unique(vapply(fits, function(f) f$data$N, 1))) > 1) cat("Note: different numbers of persons.\n")
    print_table(tab, digits)
    cat("BF_best: Bayes factor of the best model against each model (marginal likelihoods depend on the priors).\n")
    return(invisible(tab))
  }
  tab <- do.call(rbind, lapply(names(fits), function(m) cbind(model = m, fits[[m]]$fit_indices)))
  tab$dAIC <- tab$AIC - min(tab$AIC); tab$dBIC <- tab$BIC - min(tab$BIC)
  if (length(unique(vapply(fits, function(f) f$data$N, 1))) > 1) cat("Note: different numbers of persons.\n")
  print_table(tab[, c("model", "method", "npar", "Deviance", "AIC", "BIC", "dAIC", "dBIC")], digits)
  invisible(tab)
}

#' Plots of a JAGS fit
#'
#' `type = "loadings"`: posterior means and 95% intervals of the cross-loadings (or, without cross-loadings, of
#' all free loadings), with spike-and-slab inclusion probabilities when present. `"trace"`:
#' trace plots of the parameters with the largest R-hat (or those matching `pars`).
#' @param x A `birt` object.
#' @param type `"loadings"`, `"trace"`, or for fits with ordinal indicators `"dif"` (expected
#'   item scores at low and high values of the first moderator, for the flagged items or
#'   `items`) and `"information"` (test information, see [ordinal_information()]).
#' @param items Items for `type = "dif"`.
#' @param pars Regular expression of JAGS node names (`"trace"`).
#' @param ... Passed to [graphics::plot()].
#' @export
plot.birt <- function(x, type = c("loadings", "trace", "dif", "information"), pars = NULL, items = NULL, ...) {
  type <- match.arg(type)
  if (type %in% c("dif", "information")) return(plot_ordinal(x, type, items, ...))
  if (type == "loadings") {
    e <- x$estimates; ld <- e[e$op == "=~" & !e$fixed, ]
    first <- x$spec$ld$first[match(paste(ld$lhs, ld$rhs), paste(x$spec$ld$lhs, x$spec$ld$rhs))]
    if (any(!first)) ld <- ld[!first, ]
    K <- nrow(ld); lab <- paste(ld$lhs, "=~", ld$rhs)
    op <- graphics::par(mar = c(4, max(6, 0.55 * max(nchar(lab))), 2, if (any(!is.na(ld$p_incl))) 4 else 1)); on.exit(graphics::par(op))
    graphics::plot(ld$est, K:1, xlim = range(ld$q025, ld$q975, 0), yaxt = "n", pch = 19, xlab = "loading", ylab = "",
                   main = if (any(!first)) "Cross-loadings" else "Loadings", ...)
    graphics::segments(ld$q025, K:1, ld$q975, K:1); graphics::abline(v = 0, lty = 3)
    graphics::axis(2, K:1, lab, las = 1, cex.axis = 0.8)
    if (any(!is.na(ld$p_incl))) graphics::axis(4, K:1, sprintf("%.2f", ld$p_incl), las = 1, cex.axis = 0.7, tick = FALSE)
    return(invisible(ld))
  }
  sn <- x$summary_nodes
  sel <- if (is.null(pars)) utils::head(sn$node[order(-sn$rhat)], 9) else utils::head(sn$node[grepl(pars, sn$node)], 12)
  if (!length(sel)) stop("no parameters match '", pars, "'")
  nc <- ceiling(sqrt(length(sel)))
  op <- graphics::par(mfrow = c(ceiling(length(sel) / nc), nc), mar = c(3, 3, 2, 1)); on.exit(graphics::par(op))
  for (n in sel) {
    D <- sapply(x$draws, function(m) m[, n])
    graphics::matplot(D, type = "l", lty = 1, col = seq_len(ncol(D)), xlab = "", ylab = "",
                      main = sprintf("%s (Rhat %.2f)", n, sn$rhat[sn$node == n]), ...)
  }
  invisible(sel)
}

#' Keep a subset of the MCMC chains
#'
#' Recomputes the parameter table, scores, reliability and item tables from the kept
#' chains, e.g. after a chain was found in a separate posterior mode (with free
#' cross-loadings, ability can rotate into a second speed factor in one chain). The
#' marginal-likelihood criteria are dropped; recompute them with [add_fit_indices()].
#' @param object A `birt` fit.
#' @param chains Indices of the chains to keep.
#' @return The fit with the reduced draws.
#' @examples
#' \donttest{
#' if (requireNamespace("rjags", quietly = TRUE)) {
#'   fit <- birt("ability =~ y1 + y2 + y3 + y4 + y5", sim_rtirt(N = 200, K = 5),
#'               engine = "jags", n_chains = 3, n_iter = 1500, n_burn = 500, fit_indices = FALSE,
#'               progress = FALSE)
#'   # e.g. chain 3 in another mode: compare per-chain means, then drop it
#'   sapply(1:3, function(ch) mean(fit$draws[[ch]][, 1]))
#'   fit <- subset_chains(fit, c(1, 2))
#' }
#' }
#' @export
subset_chains <- function(object, chains) {
  if (!inherits(object, "birt_fit")) stop("subset_chains() takes a birt or hgrm fit")
  nc <- coda::nchain(object$draws)
  if (!length(chains) || any(!chains %in% seq_len(nc))) stop("chains must be in 1..", nc)
  object$draws <- object$draws[chains]
  object$mcmc$n_chains <- length(chains)
  object$mcmc$kept_chains <- chains
  object[c("fit_indices", "loglik", "loo", "waic")] <- list(NULL)
  add_derived(object)
}
