#' Moderated (heteroscedastic) graded response model
#'
#' A graded response model whose item discriminations and thresholds, and the mean and
#' variance of the latent trait, depend on person-level moderators (moderated nonlinear
#' factor analysis for ordinal items; Bauer, 2017). The typical use is one response
#' time per person (e.g. the total completion time of a questionnaire) moderating a
#' Likert scale: slower or faster respondents may answer with different precision
#' (discrimination, Ferrando, 2009) or shifted thresholds, a DIF-like effect of time.
#' `hgrm()` is a shortcut: it writes the model in the syntax of [birt()] with
#' `ordered = TRUE` (stored in `fit$model`), where every effect can also be set item by item,
#' and fits it with `engine = "rtmb"`: maximum likelihood with per-person adaptive quadrature
#' (marginal likelihood), or ELGM with `control = rtmb_control(method = "elgm")`.
#'
#' \deqn{P(y_{ij} > k) = \mathrm{logit}^{-1}\{a_{ij}(\theta_j - b_{ik} - \delta_{ij})\}}
#' \deqn{\log a_{ij} = \log a_i + \sum_m \alpha_{im} z_{jm}, \quad
#'       \delta_{ij} = \sum_m \beta_{im} z_{jm}}
#' \deqn{\theta_j \sim N\{\textstyle\sum_m \gamma_m z_{jm},\ \exp(2 \sum_m \psi_m z_{jm})\}}
#'
#' @param model One-factor syntax, `"f =~ y1 + y2 + ..."`. Items are ordinal with
#'   integer categories; each item's observed categories are recoded to 1..K.
#' @param data A data frame.
#' @param moderators Names of numeric person-level moderators (e.g. `"logT"`); several
#'   are allowed (e.g. a log time and its square, or piecewise terms). Standardized
#'   before fitting unless `standardize = FALSE`, so effects are per SD.
#' @param a,b Moderation of discriminations (`a`: alpha) and thresholds (`b`: beta):
#'   `"free"` one effect per item (DIF), `"common"` one effect shared by all items (for `a`:
#'   a person-level precision effect as in Ferrando, 2009), or `"none"`. A vector of length
#'   `length(moderators)` sets one option per moderator.
#' @param impact Which latent moments depend on the moderators: any of `"mean"`,
#'   `"var"`, or `"none"`.
#' @param anchor DIF-free anchor items: no threshold shifts and no item-specific (`"free"`)
#'   discrimination effects. A common discrimination effect (`a = "common"`) is a person-level
#'   precision effect and applies to every item.
#' @param standardize Standardize the moderators.
#' @param control RTMB settings, [rtmb_control()]; `method = "elgm"` for the posterior.
#' @param dp Priors of the ELGM fit, [dpriors()], as in [birt()].
#' @param ... Passed to [birt()] (`progress`, `label`, `itemtype`, ...).
#' @details
#' **Identification.** At z = 0 the trait is N(0, 1). A common threshold shift of all
#' items is the same as a shift of the trait mean, so `impact = "mean"` cannot be combined
#' with `b = "common"` or `b = "free"` unless anchors are given. Likewise a common change of
#' all discriminations and the trait variance differ only through the scaling of the
#' thresholds, so `a = "common"` or `"free"` together with `impact = "var"` is only weakly
#' identified.
#'
#' **DIF.** Fit the model without threshold effects and test them with [score_test()]
#' (one fit, every item), then free the flagged items' shifts with the others as anchors:
#' `b = "free", anchor = ...`. Spike-and-slab screening (`prior("ssp")`) is written in the
#' syntax of [birt()] with `engine = "jags"`.
#'
#' **Priors** (ELGM only; maximum likelihood has none): those of [birt()] ([dpriors()]):
#' discriminations hierarchical lognormal, first threshold `N(0, 2^2)`, increments lognormal
#' `(-0.5, 1)`, `gamma ~ N(0, 1)`, `psi ~ N(0, 1)`.
#' @return A [birt()] fit (`fit$model` holds the syntax). Fits with moderation also have
#'   `dif` (one row per moderated item parameter), `impact` and `precision` (by moderator
#'   tercile); see also [ordinal_information()] and `plot(fit, type = "dif")`.
#' @examples
#' \donttest{
#' if (requireNamespace("RTMB", quietly = TRUE)) {
#' d <- sim_hgrm(N = 400, J = 6)
#' fit <- hgrm("theta =~ y1 + y2 + y3 + y4 + y5 + y6", data = d, moderators = "logT",
#'             a = "common", impact = "mean", control = rtmb_control(method = "aghq"),
#'             progress = FALSE)
#' summary(fit)
#' fit$model                      # the same model in birt() syntax
#' score_test(fit, d$logT, by_item = TRUE)  # DIF along the time, item by item
#' }
#' }
#' @export
hgrm <- function(model, data, moderators, a = "common", b = "none", impact = "mean",
                 anchor = NULL, standardize = TRUE, control = rtmb_control(), dp = dpriors(), ...) {
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  miss <- setdiff(moderators, names(data)); if (length(miss)) stop("moderators not in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  if (standardize) for (m in moderators) {
    if (!is.numeric(data[[m]]) || anyNA(data[[m]])) stop("moderators must be numeric without missing values", call. = FALSE)
    if (stats::sd(data[[m]]) == 0) stop("constant moderator: ", m, call. = FALSE)
    data[[m]] <- as.numeric(scale(data[[m]]))
  }
  fit <- birt(hgrm_syntax(model, moderators, a, b, impact, anchor), data, ordered = TRUE, engine = "rtmb", control = control, dp = dp, ...)
  fit$call <- match.call()
  fit
}

# hgrm() options -> birt syntax: E(items) ~ z (thresholds), E(items) ~ z:f (discriminations), E(f) ~ z, V(f) ~ z
hgrm_syntax <- function(model, moderators, a, b, impact, anchor) {
  pt <- lavaan::lavParseModelString(model, as.data.frame. = TRUE, parser = "old")
  ld <- pt[pt$op == "=~", ]
  if (!nrow(ld) || length(unique(ld$lhs)) != 1 || nrow(pt) != nrow(ld))
    stop("hgrm() takes a one-factor model 'f =~ y1 + y2 + ...' (no other operators)")
  if (any(ld$fixed != "" | ld$prior != "" | ld$label != "")) stop("hgrm(): modifiers are not supported")
  items <- ld$rhs; fac <- ld$lhs[1]; M <- length(moderators)
  if (!M) stop("give at least one moderator")
  opts <- c("free", "common", "none")
  a <- rep_len(a, M); b <- rep_len(b, M)
  if (!all(c(a, b) %in% opts)) stop("a and b must be one of: ", paste(opts, collapse = ", "))
  if (length(setdiff(impact, c("mean", "var", "none")))) stop("impact must be 'mean', 'var' (or both) or 'none'")
  if (!is.null(anchor) && length(setdiff(anchor, items))) stop("anchor items not in the model: ", paste(setdiff(anchor, items), collapse = ", "))
  free <- setdiff(items, anchor)
  lines <- sprintf("%s =~ %s", fac, paste(items, collapse = " + "))
  if (!length(free) && any(c(a[a != "common"], b) != "none")) stop("hgrm(): every item is an anchor; leave at least one item free, or a = b = \"none\"", call. = FALSE)
  for (m in seq_len(M)) for (par in c("b", "a")) {
    o <- if (par == "a") a[m] else b[m]
    if (o == "none") next
    mod <- switch(o, free = "", common = sprintf("c%s_%s*", par, safe(moderators[m])))
    tg <- if (par == "a" && o == "common") items else free           # a common discrimination effect is person-level: every item
    lines <- c(lines, sprintf("E(%s) ~ %s%s%s", paste(tg, collapse = " + "), mod, moderators[m], if (par == "a") paste0(":", fac) else ""))
  }
  if ("mean" %in% impact) lines <- c(lines, sprintf("E(%s) ~ %s", fac, paste(moderators, collapse = " + ")))
  if ("var" %in% impact) lines <- c(lines, sprintf("V(%s) ~ %s", fac, paste(moderators, collapse = " + ")))
  paste(lines, collapse = "\n")
}

# ---- outputs of fits with ordinal indicators or moderation -----------------------------------

# moderator values (N x Mz) of a fit
mod_values <- function(fit) {
  zn <- fit$spec$modvars %||% character()
  if (!length(zn)) return(matrix(0, fit$data$N %||% fit$S$N, 0))
  if (inherits(fit, "birt_rtmb")) fit$S$X[, zn, drop = FALSE]
  else matrix(unlist(fit$data[zn]), fit$data$N, length(zn), dimnames = list(NULL, zn))
}

# moderation table (one row per moderated loading or shift), impact table (latent mean and log SD on
# the moderators) and the precision of the scores by tercile of the first moderator
mod_outputs <- function(fit) {
  e <- fit$estimates; sp <- fit$spec; ml <- inherits(fit, "birt_rtmb") && is.null(fit$rtmb$elgm)
  col <- function(n, alt) if (n %in% names(e)) e[[n]] else if (alt %in% names(e)) e[[alt]] else rep(NA_real_, nrow(e))
  sd <- col("sd", "se"); lo <- col("q025", "ci.lower"); hi <- col("q975", "ci.upper")
  i <- which(e$op == "~" & grepl("^E\\(", e$lhs) & !sub("^E\\((.*)\\)$", "\\1", e$lhs) %in% sp$factors)
  if (length(i)) {
    it <- sub("^E\\((.*)\\)$", "\\1", e$lhs[i]); al <- grepl(":", e$rhs[i])
    mo <- sp$mod
    ty <- vapply(seq_along(i), function(k) { r <- mo[mo$target == it[k] & mo$kind == (if (al[k]) "alpha" else "beta") & mo$mod == sub(":.*", "", e$rhs[i[k]]), ]
      if (!nrow(r)) "free" else switch(r$type[1], none = "anchor", ssp = "ssp", common = sprintf("common (%s)", r$label[1]), free = if (nzchar(r$prior[1])) r$prior[1] else "free") }, "")
    fit$dif <- data.frame(item = it, parameter = ifelse(al, "loading (log scale)", "shift"), moderator = sub(":.*", "", e$rhs[i]),
                          prior = ty, est = e$est[i], sd = sd[i], q025 = lo[i], q975 = hi[i], stringsAsFactors = FALSE)
    if (ml) { fit$dif$z <- e$z[i]; fit$dif$pvalue <- e$pvalue[i] }
    else { fit$dif$p_incl <- col("p_incl", "p_incl")[i]; fit$dif$BF10 <- col("BF10", "BF10")[i] }
    if (ml) names(fit$dif)[names(fit$dif) == "sd"] <- "se"
  }
  zn <- sp$modvars %||% character()
  j <- which((e$op == "~" & e$lhs %in% sp$factors & e$rhs %in% zn) | (e$op == "~" & grepl("^V\\(", e$lhs) & sub("^V\\((.*)\\)$", "\\1", e$lhs) %in% sp$factors))
  if (length(j)) fit$impact <- data.frame(latent = sub("^V\\((.*)\\)$", "\\1", e$lhs[j]), parameter = ifelse(grepl("^V\\(", e$lhs[j]), "log SD", "mean"),
                                          moderator = e$rhs[j], est = e$est[j], sd = sd[j], q025 = lo[j], q975 = hi[j], stringsAsFactors = FALSE)
  Z <- mod_values(fit)
  if (ncol(Z)) {
    f <- sp$factors[1]; sc <- fit$scores
    z <- Z[, 1]; br <- unique(stats::quantile(z, seq(0, 1, length.out = 4))); g <- cut(z, br, include.lowest = TRUE)
    fit$precision <- do.call(rbind, lapply(levels(g), function(l) {
      k <- which(g == l); m <- sc[[f]][k]; s <- sc[[paste0(f, "_psd")]][k]
      data.frame(group = l, n = length(k), mean_psd = mean(s), reliability = var_rel(m, s), row.names = NULL)
    }))
    names(fit$precision)[1] <- colnames(Z)[1]
  }
  fit
}

# "logT: shift ssp 8, anchor 2; loading common 10" for the header
moderation_txt <- function(sp) {
  mo <- sp$mod; if (is.null(mo) || !nrow(mo)) return("latent mean / SD only")
  out <- character()
  for (z in unique(mo$mod)) for (kd in c("beta", "alpha")) {
    r <- mo[mo$mod == z & mo$kind == kd, ]; if (!nrow(r)) next
    tb <- table(factor(r$type, c("ssp", "free", "common", "none"), c("ssp", "free", "equal", "anchor")))
    out <- c(out, sprintf("%s %s: %s", if (kd == "beta") "shift" else "loading", z, paste(sprintf("%s %d", names(tb)[tb > 0], tb[tb > 0]), collapse = ", ")))
  }
  paste(out, collapse = "; ")
}

# summary() table of the moderation effects (flagged: BF10 >= 3 or an interval excluding 0; ML: p < .05)
mod_summary <- function(fit) {
  d <- fit$dif; if (is.null(d)) return(NULL)
  d$flagged <- ifelse(rownames(d) %in% rownames(dif_flags(d)), "*", "")
  if ("p_incl" %in% names(d) && all(is.na(d$p_incl))) d$p_incl <- d$BF10 <- NULL
  btable(d, "moderation effects of the indicators; flagged: BF10 >= 3, or 95% interval excluding 0 (ML: p < .05)")
}

# flagged moderation effects: BF10 >= 3, or a 95% interval excluding 0 (Bayes); p < .05 (ML)
dif_flags <- function(dif) {
  if (is.null(dif)) return(data.frame())
  keep <- if ("pvalue" %in% names(dif)) !is.na(dif$pvalue) & dif$pvalue < .05
          else (!is.na(dif$BF10) & dif$BF10 >= 3) | (is.na(dif$p_incl) & (dif$q025 > 0 | dif$q975 < 0))
  dif[keep & dif$prior != "anchor" & !startsWith(dif$prior, "common"), ]
}

# point estimates of the ordinal items: loading a, thresholds b, log-loading (alpha) and shift (beta)
# effects per moderator
ord_point <- function(fit) {
  e <- fit$estimates; sp <- fit$spec; O <- sp$ordinal %||% character(); zn <- sp$modvars %||% character()
  get <- function(l, o, r) { i <- which(e$lhs == l & e$op == o & e$rhs == r); if (length(i)) e$est[i[1]] else 0 }
  fac <- vapply(O, function(v) sp$ld$lhs[sp$ld$rhs == v & !sp$ld$zero][1], "")
  a <- vapply(O, function(v) get(fac[[v]], "=~", v), 0)
  list(items = O, factor = fac, a = a, type = sp$ord$type[O],
       r = lapply(stats::setNames(O, O), function(v) { K <- sp$ord$K[match(v, sp$ord$items)]
         c(1, vapply(seq_len(K - 2) + 1, function(k) { i <- which(e$lhs == v & e$op == "|a" & e$rhs == paste0("a", k)); if (length(i)) e$est[i] / a[[v]] else 1 }, 0)) }),
       b = lapply(stats::setNames(O, O), function(v) vapply(seq_len(sp$ord$K[match(v, sp$ord$items)] - 1), function(k) get(v, "|", paste0("t", k)), 0)),
       alpha = matrix(vapply(zn, function(z) vapply(O, function(v) get(sprintf("E(%s)", v), "~", paste0(z, ":", fac[[v]])), 0), numeric(length(O))), length(O)),
       beta = matrix(vapply(zn, function(z) vapply(O, function(v) get(sprintf("E(%s)", v), "~", z), 0), numeric(length(O))), length(O)),
       moderators = zn, Z = mod_values(fit), levels = sp$ord$levels[O])
}

#' Test information of the ordinal indicators
#'
#' For fits with ordinal indicators (graded response or partial credit), the test information of those
#' indicators as a function of the latent variable, at chosen values of the moderators.
#' @param object A [birt()] fit with ordinal indicators (either engine).
#' @param at Named list of moderator values (default: the 10th, 50th and 90th percentile of
#'   the first moderator, the others at their mean).
#' @param theta Grid of values of the latent variable.
#' @return Data frame with `theta`, the moderator values and the test information
#'   (point estimates of the item parameters).
#' @export
ordinal_information <- function(object, at = NULL, theta = seq(-3, 3, by = 0.1)) {
  P <- ord_point(object)
  if (!length(P$items)) stop("ordinal_information() is for fits with ordinal indicators", call. = FALSE)
  M <- length(P$moderators)
  if (!M) at <- list()
  else if (is.null(at)) at <- stats::setNames(list(unname(stats::quantile(P$Z[, 1], c(.1, .5, .9)))), P$moderators[1])
  grid <- if (length(at)) expand.grid(at) else data.frame(row.names = 1)
  for (m in setdiff(P$moderators, names(grid))) grid[[m]] <- mean(P$Z[, m])
  out <- NULL
  for (r in seq_len(nrow(grid))) {
    z <- if (M) unlist(grid[r, P$moderators]) else numeric(0)
    info <- rep(0, length(theta))
    for (j in seq_along(P$items)) {
      aj <- P$a[j] * exp(sum(P$alpha[j, ] * z)); dj <- sum(P$beta[j, ] * z)
      info <- info + cat_info(P$type[[j]], aj, theta - dj, P$b[[j]], P$r[[j]])
    }
    out <- rbind(out, if (ncol(grid)) data.frame(theta = theta, grid[rep(r, length(theta)), , drop = FALSE], information = info, row.names = NULL)
                      else data.frame(theta = theta, information = info))
  }
  out
}

# plot(fit, type = "information" / "dif") for fits with ordinal indicators
plot_ordinal <- function(x, type, items = NULL, ...) {
  P <- ord_point(x)
  if (!length(P$items)) stop(sprintf("type = \"%s\" is for fits with ordinal indicators", type), call. = FALSE)
  if (type == "information") {
    inf <- ordinal_information(x)
    if (!length(P$moderators)) { graphics::plot(inf$theta, inf$information, type = "l", lwd = 2, xlab = P$factor[1], ylab = "test information", main = "Test information", ...)
                                 return(invisible(inf)) }
    m <- P$moderators[1]; lv <- unique(inf[[m]])
    graphics::plot(inf$theta, inf$information, type = "n", xlab = P$factor[1], ylab = "test information", main = sprintf("Test information by %s", m), ...)
    for (k in seq_along(lv)) graphics::lines(inf$theta[inf[[m]] == lv[k]], inf$information[inf[[m]] == lv[k]], lty = k, lwd = 2)
    graphics::legend("topright", legend = sprintf("%s = %.2f", m, lv), lty = seq_along(lv), lwd = 2, bty = "n")
    return(invisible(inf))
  }
  fl <- dif_flags(x$dif); items <- intersect(items %||% unique(fl$item), P$items)
  if (!length(items)) { message("no ordinal items with flagged moderation to plot"); return(invisible(NULL)) }
  zq <- stats::quantile(P$Z[, 1], c(.1, .9)); th <- seq(-3, 3, by = 0.05)
  nc <- ceiling(sqrt(length(items)))
  op <- graphics::par(mfrow = c(ceiling(length(items) / nc), nc)); on.exit(graphics::par(op))
  for (it in items) {
    j <- match(it, P$items); lv <- P$levels[[it]]
    graphics::plot(range(th), range(lv), type = "n", xlab = P$factor[j], ylab = "expected score", main = it, ...)
    for (k in 1:2) {
      z <- colMeans(P$Z); z[1] <- zq[k]
      aj <- P$a[j] * exp(sum(P$alpha[j, ] * z)); dj <- sum(P$beta[j, ] * z)
      graphics::lines(th, drop(cat_probs(P$type[[j]], aj, th - dj, P$b[[j]], P$r[[j]]) %*% lv), lty = k, lwd = 2)
    }
    graphics::legend("topleft", legend = sprintf("%s = %.2f", P$moderators[1], zq), lty = 1:2, lwd = 2, bty = "n", cex = 0.8)
  }
  invisible(items)
}

# posterior predictive check of the ordinal indicators: category proportions by moderator tercile
ppc_ordinal <- function(P, sp) {
  O <- sp$ordinal; N <- P$N
  th_cols <- sprintf("%s[%d]", P$f1, seq_len(N))
  g <- if (ncol(P$Z)) cut(P$Z[, 1], unique(stats::quantile(P$Z[, 1], 0:3 / 3)), include.lowest = TRUE) else factor(rep("all", N))
  rows <- list()
  for (v in O) {
    K <- sp$ord$K[match(v, sp$ord$items)]; y <- sp$ord$Y[, v]
    obs <- vapply(levels(g), function(l) tabulate(y[g == l], K) / sum(!is.na(y[g == l])), numeric(K))
    rep <- array(NA, c(P$S, K, nlevels(g)))
    for (s in seq_len(P$S)) {
      fv <- sp$ld$lhs[sp$ld$rhs == v & !sp$ld$zero][1]                   # the latent variable of the item
      th <- P$M[s, sprintf("%s[%d]", fv, seq_len(N))]; a <- (P$o1[s, v] + P$o2[s, v]) * mod_mult(P, s, paste(fv, v)); d <- mod_shift(P, s, v)
      Pc <- cat_probs(P$otype[[v]], a, th - d, P$thr[[v]][s, ], P$rat[[v]][s, ])
      yr <- 1 + rowSums(stats::runif(N) > t(apply(Pc, 1, cumsum))[, -K, drop = FALSE]); yr[is.na(y)] <- NA
      rep[s, , ] <- vapply(levels(g), function(l) tabulate(yr[g == l], K) / sum(!is.na(yr[g == l])), numeric(K))
    }
    for (gi in seq_len(nlevels(g))) for (k in seq_len(K)) {
      pv <- mean(rep[, k, gi] >= obs[k, gi])
      rows[[length(rows) + 1]] <- data.frame(item = v, group = levels(g)[gi], category = sp$ord$levels[[v]][k],
                                             p_obs = obs[k, gi], p_rep = mean(rep[, k, gi]), PPP = pv, flag = ifelse(pv < .025 | pv > .975, "!", ""))
    }
  }
  tab <- do.call(rbind, rows); names(tab)[2] <- if (ncol(P$Z)) colnames(P$Z)[1] else "group"
  tab
}

#' Simulate Likert responses moderated by a total response time
#'
#' GRM with `J` items and 4 categories; the log total time `logT` lowers the
#' discrimination of every item (person precision, `alpha`), shifts the thresholds of
#' items 1 and 2 (`beta`), and is related to the trait mean (`gamma`).
#' @param N,J Persons, items.
#' @param alpha Common effect of standardized `logT` on log discrimination.
#' @param beta Threshold shift of items 1 and 2 per SD of `logT`.
#' @param gamma Effect of `logT` on the trait mean.
#' @param seed Random seed.
#' @return A data frame with `y1..yJ` (1..4), `Time` and `logT`.
#' @export
sim_hgrm <- function(N = 800, J = 10, alpha = 0.25, beta = 0.4, gamma = 0.2, seed = 1) {
  local_seed(seed)
  logT <- stats::rnorm(N, 6.3, 0.7); z <- (logT - mean(logT)) / stats::sd(logT)
  theta <- gamma * z + stats::rnorm(N)
  a <- stats::runif(J, 1, 2); b <- t(apply(matrix(stats::rnorm(J * 3, 0, 1), J), 1, sort)) + c(-0.5, 0, 0.5)
  Y <- matrix(NA, N, J)
  for (j in 1:J) {
    aij <- a[j] * exp(alpha * z); dj <- if (j <= 2) beta * z else 0
    Ps <- stats::plogis(aij * (theta - dj) - outer(aij, b[j, ]))
    Y[, j] <- 1 + rowSums(stats::runif(N) < Ps)
  }
  colnames(Y) <- paste0("y", 1:J)
  data.frame(Y, Time = round(exp(logT)), logT = logT)
}
