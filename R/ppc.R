#' Posterior predictive checks
#'
#' Replicates the data from posterior draws (person parameters included) and compares:
#' item proportions correct, the sum-score distribution, means and SDs of the continuous
#' indicators, and the conditional dependence between each binary item and its paired
#' continuous indicator (correlation of accuracy with the standardized residual;
#' Bolsinova & Tijmstra, 2016). A continuous indicator is paired with a binary item when
#' its name ends with the item name (e.g. `lrt_ME1` and `ME1`), or by position when the
#' numbers of binary and continuous indicators are equal.
#' @param object A `birt` object.
#' @param n_draws Number of posterior draws.
#' @param seed Random seed.
#' @param ... Unused.
#' @return An object of class `birt_ppc` (list of tables); PPP values outside
#'   (.025, .975) are flagged.
#' @export
ppc <- function(object, ...) UseMethod("ppc")

#' @rdname ppc
#' @export
ppc.birt <- function(object, n_draws = 500, seed = 1, ...) {
  local_seed(seed)
  P <- prep_pars(object, n_draws); sp <- object$spec; N <- P$N; M <- P$M
  th <- M[, sprintf("%s[%d]", P$f1, seq_len(N)), drop = FALSE]
  s2 <- if (!is.na(P$f2)) M[, sprintf("%s[%d]", P$f2, seq_len(N)), drop = FALSE] else matrix(0, P$S, N)
  ob <- obs_matrices(object); Y <- ob$Y; Tm <- ob$T
  Kb <- ncol(Y); C <- ncol(Tm); full <- Kb > 0 && !anyNA(Y)
  pr <- pairs_of(sp)
  cor_na <- function(a, b) vapply(seq_len(ncol(a)), function(k) suppressWarnings(stats::cor(a[, k], b[, k], use = "complete.obs")), 0)
  r_p <- matrix(NA, P$S, Kb); r_m <- r_s <- matrix(NA, P$S, C)
  o_cd <- r_cd <- matrix(NA, P$S, length(pr)); r_score <- matrix(NA, P$S, Kb + 1)
  for (s in seq_len(P$S)) {
    yr <- NULL
    if (Kb) {
      eta <- vapply(sp$binary, function(v) {                       # loadings on both latent variables, with moderation
        l1 <- P$a1[s, v] * mod_mult(P, s, paste(P$f1, v)); l2 <- if (is.na(P$f2)) 0 else P$a2[s, v] * mod_mult(P, s, paste(P$f2, v))
        lf <- if (P$a1[s, v] != 0) l1 else l2                          # the shift is on the scale of the (first) loading
        rep_len(P$d[s, v] + l1 * th[s, ] + l2 * s2[s, ] - lf * mod_shift(P, s, v), N)
      }, numeric(N))
      eta <- matrix(eta, N)
      yr <- matrix(stats::rbinom(N * Kb, 1, stats::plogis(eta)), N, Kb, dimnames = list(NULL, sp$binary)); yr[is.na(Y)] <- NA
      r_p[s, ] <- colMeans(yr, na.rm = TRUE)
      if (full) r_score[s, ] <- tabulate(rowSums(yr) + 1, Kb + 1)
    }
    if (C) {
      cl <- function(f) matrix(vapply(sp$cont, function(v) rep_len(f(v), N), numeric(N)), N)
      mu <- matrix(P$xi[s, ], N, C, byrow = TRUE) + cl(function(v) P$l1[s, v] * mod_mult(P, s, paste(P$f1, v))) * th[s, ] +
            cl(function(v) if (is.na(P$f2)) 0 else P$l2[s, v] * mod_mult(P, s, paste(P$f2, v))) * s2[s, ] + cl(function(v) mod_shift(P, s, v))
      sgm <- cl(function(v) P$sig[s, match(v, sp$cont)] * mod_sd(P, s, v))
      tr <- mu
      for (c in seq_len(C)) tr[, c] <- switch(P$ct[c], normal = mu[, c] + sgm[, c] * stats::rnorm(N),
                                              ald = rald(N, mu[, c], sgm[, c], P$pc[c]))
      sdc <- vapply(seq_len(C), function(c) switch(P$ct[c], normal = 1, ald = sqrt(ald_var(P$pc[c]))), 0)
      tr[is.na(Tm)] <- NA
      colnames(tr) <- colnames(mu) <- sp$cont
      r_m[s, ] <- colMeans(tr, na.rm = TRUE); r_s[s, ] <- apply(tr, 2, stats::sd, na.rm = TRUE)
      if (length(pr)) {
        sg <- sgm * matrix(sdc, N, C, byrow = TRUE); colnames(sg) <- sp$cont
        o_cd[s, ] <- cor_na(Y[, names(pr), drop = FALSE], ((Tm - mu) / sg)[, pr, drop = FALSE])
        r_cd[s, ] <- cor_na(yr[, names(pr), drop = FALSE], ((tr - mu) / sg)[, pr, drop = FALSE])
      }
    }
  }
  ppp <- function(rep, obs) colMeans(sweep(rep, 2, obs, ">="), na.rm = TRUE)
  flag <- function(p) ifelse(p < 0.025 | p > 0.975, "!", "")
  out <- list(n_draws = P$S)
  if (Kb) {
    o_p <- colMeans(Y, na.rm = TRUE); pp <- ppp(r_p, o_p)
    out$accuracy <- data.frame(item = sp$binary, p_obs = o_p, p_rep = colMeans(r_p), PPP = pp, flag = flag(pp), row.names = NULL)
  }
  if (full) {
    o_score <- tabulate(rowSums(Y) + 1, Kb + 1); E <- colMeans(r_score)
    chi <- function(o) sum((o - E)^2 / pmax(E, 0.5))
    x_r <- apply(r_score, 1, chi); pv <- mean(x_r >= chi(o_score))
    out$score <- data.frame(discrepancy = "sum-score distribution (chi-square)", obs = chi(o_score),
                            rep_mean = mean(x_r), PPP = pv, flag = flag(pv))
  }
  if (C) {
    o_m <- colMeans(Tm, na.rm = TRUE); o_s <- apply(Tm, 2, stats::sd, na.rm = TRUE)
    pm <- ppp(r_m, o_m); ps <- ppp(r_s, o_s)
    out$continuous <- data.frame(indicator = sp$cont, mean_obs = o_m, mean_rep = colMeans(r_m), PPP_mean = pm,
                                 sd_obs = o_s, sd_rep = colMeans(r_s), PPP_sd = ps, flag = paste0(flag(pm), flag(ps)),
                                 row.names = NULL)
  }
  if (length(pr)) {
    pc <- colMeans(r_cd >= o_cd, na.rm = TRUE)
    out$dependence <- data.frame(item = names(pr), indicator = unname(pr), cor_realized = colMeans(o_cd),
                                 cor_replicated = colMeans(r_cd), PPP = pc, flag = flag(pc), row.names = NULL)
  }
  if (length(sp$ordinal)) out$categories <- ppc_ordinal(P, sp)
  structure(out, class = "birt_ppc")
}

pairs_of <- function(sp) {
  if (!length(sp$cont) || !length(sp$binary)) return(character())
  m <- vapply(sp$binary, function(b) { hit <- sp$cont[endsWith(sp$cont, b)]; if (length(hit)) hit[1] else NA_character_ }, "")
  if (all(is.na(m)) && length(sp$binary) == length(sp$cont)) m <- stats::setNames(sp$cont, sp$binary)
  m[!is.na(m)]
}

#' @export
print.birt_ppc <- function(x, digits = 3, ...) {
  cat(sprintf("Posterior predictive checks (%d draws; ! = PPP < .025 or > .975)\n", x$n_draws))
  if (!is.null(x$score)) print_table(x$score, digits)
  if (!is.null(x$accuracy)) { cat("Item proportion correct\n"); print_table(x$accuracy, digits) }
  if (!is.null(x$continuous)) { cat("Continuous indicators: mean and SD\n"); print_table(x$continuous, digits) }
  if (!is.null(x$dependence)) { cat("Conditional dependence: cor(accuracy, standardized residual)\n"); print_table(x$dependence, digits) }
  if (!is.null(x$categories)) {
    fl <- x$categories[x$categories$flag != "", ]
    cat(sprintf("Category proportions by moderator tercile: %d of %d cells flagged\n", nrow(fl), nrow(x$categories)))
    if (nrow(fl)) print_table(fl, digits)
  }
  invisible(x)
}
