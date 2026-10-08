# Parameter table (lavaan-like) and item tables built from the draws ------------------

summ_vec <- function(x) c(est = mean(x), sd = stats::sd(x), q025 = unname(stats::quantile(x, 0.025)),
                          q975 = unname(stats::quantile(x, 0.975)))

# prior odds of inclusion of a spike-and-slab parameter: fixed ssp_p, or (learned p_incl ~
# beta(1, 1)) the posterior mean of the inclusion probability shared by its group (loadings
# "load", regressions "reg", moderation effects one group per kind and moderator, e.g.
# "dlt_z" / "alm_z") as a plug-in. Every ssp row of the estimates carries its group.
ssp_group <- function(r) if (!is.null(r$group) && !is.na(r$group)) r$group else if (r$kind == "loading") "load" else "reg"
prior_odds <- function(r, M, sp) {
  p <- if (!is.null(sp$ssp_p)) sp$ssp_p else {
    n <- sprintf("p_incl_%s", ssp_group(r))
    if (!n %in% colnames(M)) stop(sprintf("the inclusion probability %s was not monitored", n))
    mean(M[, n])
  }
  p / (1 - p)
}

make_estimates <- function(fit) {
  M <- as.matrix(fit$draws); pt <- fit$partab; sp <- fit$spec; sn <- fit$summary_nodes
  out <- lapply(seq_len(nrow(pt)), function(i) {
    r <- pt[i, ]
    x <- if (!is.na(r$node)) M[, r$node] else rep(r$fixed, nrow(M))
    var_kind <- grepl("^(resid|fsd|fscale)", r$kind)
    fam <- if (grepl("^resid", r$kind)) sp$family$ind[[r$lhs]]
           else if (grepl("^(fsd|fscale)", r$kind)) sp$family$factor[[r$lhs]] else NULL
    scale <- NA_real_
    if (r$kind %in% c("resid", "fsd")) x <- x^2                                   # variance
    if (r$kind %in% c("resid_ald", "fscale_ald")) {
      scale <- mean(x); x <- x^2 * ald_var(fam$p)                                 # implied variance
    }
    s <- summ_vec(x)
    fixed <- is.na(r$node) || (grepl("^fscale", r$kind) && !is.na(r$fixed))
    k <- match(r$node, sn$node)
    pin <- if (!is.na(r$incl)) mean(M[, r$incl]) else NA_real_
    data.frame(lhs = r$lhs, op = r$op, rhs = r$rhs,
               family = if (is.null(fam)) "" else fam_label(fam), prior = if (fixed) "" else r$prior,
               est = s[["est"]], sd = if (fixed) NA else s[["sd"]], q025 = if (fixed) NA else s[["q025"]],
               q975 = if (fixed) NA else s[["q975"]], rhat = if (is.na(k)) NA else sn$rhat[k],
               ess = if (is.na(k)) NA else sn$ess[k], scale = scale, p_incl = pin,
               BF10 = if (is.na(pin)) NA else (pin / (1 - pin)) / prior_odds(r, M, sp), fixed = fixed,
               kind = r$kind, node = r$node, group = if (is.na(pin)) NA_character_ else ssp_group(r), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

make_items <- function(fit) {
  P <- prep_pars(fit); sp <- fit$spec
  bin <- NULL
  if (length(sp$binary)) {
    a <- P$a1; bdraw <- -P$d / a
    bin <- data.frame(item = sp$binary, a = colMeans(a), b = colMeans(bdraw), d = colMeans(P$d),
                      row.names = NULL)
    if (any(P$a2 != 0)) bin[[sprintf("a_%s", P$f2)]] <- colMeans(P$a2)
  }
  cont <- NULL
  if (length(sp$cont)) {
    vv <- vapply(seq_along(sp$cont), function(c) switch(P$ct[c], normal = 1, ald = ald_var(P$pc[c])), 0)
    cont <- data.frame(item = sp$cont, family = vapply(sp$cont, function(v) fam_label(sp$family$ind[[v]]), ""),
                       xi = colMeans(P$xi), row.names = NULL)
    cont[[sprintf("lambda_%s", P$f1)]] <- colMeans(P$l1)
    if (!is.na(P$f2)) cont[[sprintf("lambda_%s", P$f2)]] <- colMeans(P$l2)
    cont$sigma <- colMeans(P$sig)
    cont$var <- colMeans(P$sig^2) * vv
    ic <- vapply(sp$cont, function(v) { i <- which(fit$partab$kind == "loading" & fit$partab$lhs == P$f1 & fit$partab$rhs == v)
                                        if (length(i)) fit$partab$incl[i] else NA_character_ }, "")
    if (any(!is.na(ic))) cont$p_incl <- vapply(ic, function(n) if (is.na(n)) NA_real_ else mean(P$M[, n]), 0)
  }
  ord <- NULL
  if (length(sp$ordinal)) {
    ord <- data.frame(item = sp$ordinal, K = sp$ord$K[match(sp$ordinal, sp$ord$items)], a = colMeans(P$o1 + P$o2), row.names = NULL)
    for (k in seq_len(max(ord$K) - 1)) ord[[sprintf("b%d", k)]] <- vapply(sp$ordinal, function(v) if (k < ncol(P$thr[[v]]) + 1) mean(P$thr[[v]][, k]) else NA_real_, 0)
  }
  list(binary = bin, continuous = cont, ordinal = ord)
}
