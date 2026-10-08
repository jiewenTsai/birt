# birt_rtmb objects: fit assembly and methods --------------------------------------------

birt_rtmb <- function(cl, model, data, spec, control, label, id, progress) {
  S <- rtmb_build(spec, data)
  bayes <- control$method == "elgm"
  if (bayes) S$prior <- rtmb_priors(S, spec$dp %||% dpriors())
  else {
    sp <- spec
    pri <- c(with(sp$ld, sprintf("%s =~ %s", lhs, rhs)[prior != ""]), with(sp$reg, sprintf("%s ~ %s", lhs, rhs)[prior != ""]))
    if (length(pri)) message("engine = \"rtmb\" (maximum likelihood) ignores prior() and estimates these freely: ",
                             paste(utils::head(pri, 3), collapse = ", "), if (length(pri) > 3) sprintf(", ... (%d)", length(pri)) else "",
                             "; rtmb_control(method = \"elgm\") uses them")
  }
  r <- fit_rtmb(S, control, progress)
  if (bayes) {
    lp0 <- S$prior$lp(r$par)
    r$logpost <- -r$objective; r$objective <- r$objective + lp0            # objective: -log L at the posterior mode
    r$elgm <- elgm_outer(S, r, control, progress)
    r$scores <- r$elgm$scores
  }
  fit <- structure(list(call = cl, model = model, label = label, spec = spec, S = S, rtmb = r, control = control, likelihood = "marginal",
                        id = id, data = list(N = nrow(data))), class = "birt_rtmb")
  fit$estimates <- rtmb_estimates(fit)
  e <- fit$estimates; dv <- e$op == "=~" & e$rhs %in% spec$binary & !e$fixed & abs(e$est) > 10
  if (!bayes && any(dv)) warning(sprintf("diverging discrimination of %s (loading %s): the maximum likelihood estimate does not exist for these data (the item is close to a step function of %s), and the quadrature is unreliable; fix or drop the item, or use engine = \"jags\" (the priors keep it finite)",
                               paste(e$rhs[dv], collapse = ", "), paste(format(round(e$est[dv])), collapse = ", "), e$lhs[dv][1]), call. = FALSE)
  N <- S$N; sc <- data.frame(id = id %||% seq_len(N)); rel <- numeric()
  for (i in seq_len(S$nF)) {
    f <- S$fs[i]; m <- r$scores[[i]]$mean; s <- r$scores[[i]]$sd
    sc[[f]] <- m; sc[[paste0(f, "_psd")]] <- s
    rel[f] <- var_rel(m, s)
  }
  fit$scores <- sc; fit$reliability <- rel
  fit <- mod_outputs(fit)
  ll <- -r$objective
  inner <- if (r$method == "aghq") sprintf("AGHQ (k = %d)", r$k) else "Laplace"
  fit$fit_indices <- if (bayes) {
    E <- r$elgm
    data.frame(method = sprintf("ELGM: %s inner, %s outer", inner,
                                if (E$nq) sprintf("AGHQ (k = %d) on %d", E$k, E$nq) else "Gaussian"),
               npar = r$npar, logML = E$lognc, logLik_mode = ll, Deviance_mode = -2 * ll, nodes = length(E$prob), ESS_nodes = E$ess)
  } else data.frame(method = inner, npar = r$npar, logLik = ll, Deviance = -2 * ll, AIC = -2 * ll + 2 * r$npar,
                    BIC = -2 * ll + log(N) * r$npar)
  if (r$convergence != 0 || r$maxgrad > 0.01)
    warning(sprintf("the optimizer may not have converged (code %d, max |gradient| %.2g); see convergence()", r$convergence, r$maxgrad), call. = FALSE)
  fit
}

#' The `birt_rtmb` object (RTMB engine)
#'
#' A model fitted with `birt(..., engine = "rtmb")`: maximum likelihood with the latent
#' variables integrated out by adaptive Gauss-Hermite quadrature (or Laplace). Accessors as
#' for JAGS fits: `print`, `summary`, `coef`, [estimates()] (est, SE, z, p, 95\% CI),
#' [scores()] (EAP and posterior SD from the quadrature), [reliability()],
#' [fit_indices()] (log-likelihood, AIC, BIC), [convergence()], `logLik`, `vcov`, `AIC`,
#' `BIC`, [quantiles()] (quantile effects) and `plot` (quantile curves).
#' @name birt_rtmb
#' @param object,x A `birt_rtmb` fit.
#' @param digits Number of decimals.
#' @param ... Unused, or passed to [graphics::plot()].
NULL

rtmb_header <- function(x) {
  sp <- x$spec; r <- x$rtmb; fi <- x$fit_indices
  c(sprintf("birt %s, RTMB engine: marginal likelihood (persons integrated out)%s", pkg_version(), if (!is.null(x$label)) sprintf(" [%s]", x$label) else ""),
    model_lines(x),
    if (length(sp$log_rt)) sprintf("  Log-transformed: %s", paste(c(utils::head(sp$log_rt, 4), if (length(sp$log_rt) > 4) sprintf("... (%d)", length(sp$log_rt))), collapse = ", ")),
    censor_txt(sp$censor),
    if (is.null(r$elgm)) c(
      sprintf("  Maximum likelihood, %s%s: %d parameters; %.0f s", fi$method,
              if (r$method == "aghq") sprintf(", %d rounds", length(r$rounds)) else "", r$npar, r$secs),
      sprintf("  -2logL %.1f, AIC %.1f, BIC %.1f", fi$Deviance, fi$AIC, fi$BIC))
    else {
      E <- r$elgm
      c(sprintf("  Approximate Bayes (ELGM; W = person latents): %d parameters; %.0f s", r$npar, r$secs + E$secs),
        sprintf("    inner: %s per person; outer: %s", if (E$inner == "aghq") sprintf("AGHQ (k = %d)", E$k_inner) else "Laplace",
                if (E$nq) sprintf("AGHQ, k = %d on %s, Gaussian in %d directions (%d nodes)", E$k, paste(E$hyper, collapse = ", "), E$n_gauss, length(E$prob))
                else "Gaussian at the posterior mode (Laplace)"),
        sprintf("  log marginal likelihood %.2f; -2logL at the posterior mode %.1f", E$lognc, fi$Deviance_mode))
    })
}

#' @rdname birt_rtmb
#' @export
print.birt_rtmb <- function(x, ...) {
  cat(rtmb_header(x), sep = "\n")
  cat(sprintf("  Reliability: %s\n", paste(sprintf("%s %.3f", names(x$reliability), x$reliability), collapse = ", ")))
  if (x$rtmb$convergence != 0 || x$rtmb$maxgrad > 0.01) cat("  Warning: check convergence()\n")
  cat("Use summary() for the parameter table and quantiles() / plot() for the quantile effects.\n")
  invisible(x)
}

#' @rdname birt_rtmb
#' @export
summary.birt_rtmb <- function(object, digits = 3, ...) {
  fit <- object; e <- fit$estimates; sp <- fit$spec; bayes <- !is.null(fit$rtmb$elgm)
  if (bayes) { names(e)[match(c("se", "ci.lower", "ci.upper"), names(e))] <- c("sd", "q025", "q975") }
  par <- partable(e[!e$op %in% c("irt_a", "irt_b"), ], bayes,
                  if (bayes) "approximate posterior (ELGM): mean, SD and 95% interval; prior of each parameter"
                  else "maximum likelihood: estimate, standard error, Wald z and 95% interval",
                  sections = c("Latent Variables", "Regressions", "Moderation (E: expectation, V: variance)", "Covariances (correlations)", "Variances",
                               "Distributions", "Intercepts", "Thresholds"), spec = sp)
  pe <- as.data.frame(unclass(partable(e, bayes, "")))[, c("lhs", "op", "rhs", "est", "se")]
  w <- item_wide(pe, sp$ovs, sp$factors, sp$binary, sp$ordinal)
  structure(list(header = rtmb_header(fit),
                 fit = btable(cbind(N = fit$data$N, fit$fit_indices), if (bayes) "log marginal likelihood (ELGM)" else "log-likelihood (persons integrated out by AGHQ)"),
                 items = btable(w$est, sprintf("item parameters (%s)%s", if (bayes) "posterior means" else "ML", if (length(sp$binary)) "; binary items also as IRT a, b = -d / a" else "")),
                 se = btable(w$se, if (bayes) "posterior SDs of the item parameters" else "standard errors of the item parameters"),
                 parameters = par, moderation = mod_summary(fit),
                 precision = btable(fit$precision, "precision of the scores by tercile of the first moderator"),
                 reliability = btable(as.data.frame(as.list(fit$reliability)), "reliability of the EAP scores: var(EAP) / (var(EAP) + mean PSD^2)")),
            class = c("summary.birt_rtmb", "birt_summary"))
}

#' @rdname birt_rtmb
#' @export
estimates.birt_rtmb <- function(object, ...) {
  e <- object$estimates; e <- e[!e$op %in% c("irt_a", "irt_b"), ]
  if (is.null(object$rtmb$elgm)) return(e[, c("lhs", "op", "rhs", "est", "se", "z", "pvalue", "ci.lower", "ci.upper")])
  # approximate Bayes: the column names of JAGS fits (posterior mean, SD, 2.5% and 97.5%)
  out <- e[, c("lhs", "op", "rhs", "prior", "est", "se", "ci.lower", "ci.upper")]
  names(out)[6:8] <- c("sd", "q025", "q975")
  out
}

#' @rdname birt_rtmb
#' @export
coef.birt_rtmb <- function(object, ...) {
  e <- object$estimates; e <- e[!e$fixed & !e$op %in% c("irt_a", "irt_b"), ]
  stats::setNames(e$est, ifelse(e$op %in% c("=~", "~", "~~", "~1"), paste0(e$lhs, e$op, e$rhs), paste0(e$lhs, ":", e$op, e$rhs)))
}

#' @rdname birt_rtmb
#' @export
vcov.birt_rtmb <- function(object, ...) {
  V <- object$rtmb$vcov; dimnames(V) <- list(names(object$rtmb$x), names(object$rtmb$x)); V
}

#' @rdname birt_rtmb
#' @export
logLik.birt_rtmb <- function(object, ...)
  structure(-object$rtmb$objective, df = object$rtmb$npar, nobs = object$data$N, class = "logLik")

#' @rdname birt_rtmb
#' @export
scores.birt_rtmb <- function(object, ...) object$scores
#' @rdname birt_rtmb
#' @export
reliability.birt_rtmb <- function(object, ...) object$reliability
#' @rdname birt_rtmb
#' @export
fit_indices.birt_rtmb <- function(object, ...) object$fit_indices

#' @rdname birt_rtmb
#' @export
convergence.birt_rtmb <- function(object, ...) {
  r <- object$rtmb
  tab <- data.frame(method = object$fit_indices$method, nlminb_code = r$convergence, max_abs_gradient = r$maxgrad,
                    se_available = mean(!is.na(object$estimates$se[!object$estimates$fixed])))
  print_table(tab, 4)
  if (!is.null(r$elgm)) {
    E <- r$elgm
    cat(sprintf("ELGM outer quadrature: %d nodes, effective %.1f; largest mass on the two outermost nodes of a direction %.3f (Gaussian: %.3f); max log posterior above the mode %.3f\n",
                length(E$prob), E$ess, E$edge, E$edge0 %||% NA, E$excess))
    cat("(the AGHQ rounds below refer to the posterior mode, i.e. -2 log L - 2 log prior)\n")
  }
  if (r$method == "aghq") {
    cat(sprintf("AGHQ rounds (-2logL): %s\n", paste(sprintf("%.3f", 2 * r$rounds), collapse = ", ")))
    cat(sprintf("Quadrature check (-2logL at the estimates): %s\n", paste(sprintf("%s %.3f", names(r$quad_check), r$quad_check), collapse = ", ")))
  }
  cat(sprintf("Laplace -2logL (start): %.3f\n", 2 * r$laplace$objective))
  invisible(tab)
}

#' @rdname birt_rtmb
#' @param type `"quantiles"`: quantile effects against the level p (one panel per target and
#'   predictor whose effect can vary with p); `"scores"`: EAP scores of the two latent variables;
#'   `"dif"`, `"information"`: ordinal indicators, as for [plot.birt()].
#' @param items Items for `type = "dif"`.
#' @param p Quantile levels.
#' @param targets,predictors Restrict the panels (names of latent variables or indicators /
#'   of predictors).
#' @param ald Optional `birt_quantile` fit (`engine = "jags"`, [ald()] at several p) whose
#'   coefficients are overlaid (points with 95\% intervals).
#' @param max_panels Maximum number of panels.
#' @export
plot.birt_rtmb <- function(x, type = c("quantiles", "scores", "dif", "information"), p = seq(0.05, 0.95, by = 0.05), targets = NULL,
                           predictors = NULL, ald = NULL, max_panels = 16, items = NULL, ...) {
  type <- match.arg(type)
  if (type %in% c("dif", "information")) return(plot_ordinal(x, type, items, ...))
  if (type == "scores") {
    sc <- x$scores; fs <- x$spec$factors
    if (length(fs) == 1) graphics::hist(sc[[fs]], main = fs, xlab = "EAP", ...)
    else graphics::plot(sc[[fs[1]]], sc[[fs[2]]], xlab = fs[1], ylab = fs[2], pch = 20, col = grDevices::adjustcolor(1, 0.4), ...)
    return(invisible(sc))
  }
  q <- quantiles(x, p)
  q <- q[!is.na(q$se), ]
  if (!nrow(q)) { message("no quantile effects to plot (no latent or continuous targets with predictors); plot(fit, type = \"scores\")"); return(invisible(q)) }
  vary <- stats::aggregate(est ~ target + predictor, q, function(v) diff(range(v)))
  key <- vary[vary$est > 1e-8, c("target", "predictor")]
  if (!is.null(targets)) key <- key[key$target %in% targets, ]
  if (!is.null(predictors)) key <- key[key$predictor %in% predictors, ]
  if (!nrow(key)) { message("no quantile effects that change with p: no predictor moderates the scale or skewness of a target, so each effect is a location shift, the same at every p; see quantiles()"); return(invisible(q)) }
  if (nrow(key) > max_panels) { message(sprintf("showing %d of %d panels; use targets = / predictors =", max_panels, nrow(key))); key <- key[seq_len(max_panels), ] }
  A <- if (!is.null(ald)) {
    if (!inherits(ald, "birt_quantile")) stop("`ald` must be a birt_quantile fit (engine = \"jags\", ald() at several p)")
    estimates(ald)
  }
  nc <- ceiling(sqrt(nrow(key)))
  op <- graphics::par(mfrow = c(ceiling(nrow(key) / nc), nc), mar = c(4, 4, 2, 1)); on.exit(graphics::par(op))
  out <- list()
  for (i in seq_len(nrow(key))) {
    tg <- key$target[i]; pr <- key$predictor[i]
    r <- q[q$target == tg & q$predictor == pr, ]; r <- r[order(r$p), ]
    a <- if (!is.null(A)) {
      if (tg %in% x$spec$factors) A[A$op == "~" & A$lhs == tg & A$rhs == pr, ] else A[A$op == "=~" & A$lhs == pr & A$rhs == tg, ]
    }
    yl <- range(r$ci.lower, r$ci.upper, if (!is.null(a) && nrow(a)) c(a$q025, a$q975), 0)
    graphics::plot(r$p, r$est, type = "n", ylim = yl, xlab = "quantile p", ylab = sprintf("effect of %s", pr), main = tg, ...)
    graphics::polygon(c(r$p, rev(r$p)), c(r$ci.lower, rev(r$ci.upper)), col = "grey85", border = NA)
    graphics::lines(r$p, r$est, lwd = 2); graphics::abline(h = 0, lty = 3)
    if (!is.null(a) && nrow(a)) {
      graphics::segments(a$p, a$q025, a$p, a$q975, col = "firebrick"); graphics::points(a$p, a$est, pch = 19, col = "firebrick")
    }
    out[[i]] <- r
  }
  invisible(do.call(rbind, out))
}
