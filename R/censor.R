# Censored continuous indicators (time limits, fast responses known only as fast) -------------
#
# A response time at or beyond its upper limit c is right censored: only "at least c" is known,
# so it contributes log S(c | f) = log P(y >= c | f) instead of the log density; a value at or
# below a lower limit l is left censored and contributes log F(l | f) (the usual likelihood of
# censored data, e.g. Meeker & Escobar, 1998; time limits in IRT: Lee & Ying, 2015). With y = mu + sd e and e = T(z), z ~ N(0, 1) (T the identity for
# normal(), the standardized SHASH transform T_std for shash()), F(c | f) = Phi(T^-1((c - mu) / sd)).

# T_std^-1: the standard normal value z with T_std(z; eps, delta) = x (see T_std() in R/rtmb.R).
# The SHASH argument w is smoothly bounded at +-40 as in ld_sh() (identity to 1e-7 for |w| < 8,
# where Phi(-sinh 8) = exp(-1.1e6) already), so log Phi stays finite at far nodes.
T_inv <- function(x, eps, dl) {
  a <- eff_a(eps, dl)
  w <- dl * (asinh(x * cosh(a) / dl + sinh(a)) - a)
  sinh(w / (1 + (w / 40)^8)^0.125)
}

# censor = list(t1 = c(upper = 60), t2 = c(lower = 1, upper = 60)): checked against the
# continuous indicators; limits are on the scale of the columns of `data` as given and are
# transformed with them (log for the columns of log_rt). Returns a data frame (one row per
# censored indicator) with the limits on the model scale and the counts, or NULL.
censor_setup <- function(censor, sp, data, logv = character()) {
  if (is.null(censor) || !length(censor)) return(NULL)
  if (!is.list(censor) || is.null(names(censor)) || any(!nzchar(names(censor))))
    stop("`censor` must be a named list, e.g. censor = list(t1 = c(upper = 60))", call. = FALSE)
  bad <- setdiff(names(censor), sp$cont)
  if (length(bad)) stop("`censor`: not continuous indicators of the model: ", paste(bad, collapse = ", "), call. = FALSE)
  if (anyDuplicated(names(censor))) stop("`censor`: an indicator is named twice", call. = FALSE)
  rows <- lapply(names(censor), function(v) {
    l <- censor[[v]]
    if (!is.numeric(l) || is.null(names(l)) || !all(names(l) %in% c("lower", "upper")) || anyDuplicated(names(l)) || anyNA(l))
      stop(sprintf("`censor$%s` must be a named number or pair: c(upper = 60), c(lower = 1), c(lower = 1, upper = 60)", v), call. = FALSE)
    lo <- if ("lower" %in% names(l)) unname(l["lower"]) else -Inf; up <- if ("upper" %in% names(l)) unname(l["upper"]) else Inf
    if (lo >= up) stop(sprintf("`censor$%s`: lower must be below upper", v), call. = FALSE)
    lg <- v %in% logv
    if (lg && any(c(lo, up)[is.finite(c(lo, up))] <= 0)) stop(sprintf("`censor$%s`: %s is log-transformed (log_rt), so its limits must be positive seconds", v, v), call. = FALSE)
    y <- data[[v]]
    LO <- if (lg && is.finite(lo)) log(lo) else lo; UP <- if (lg && is.finite(up)) log(up) else up
    data.frame(item = v, lower = LO, upper = UP, lower_raw = lo, upper_raw = up, log = lg,
               n_left = sum(!is.na(y) & y <= LO), n_right = sum(!is.na(y) & y >= UP), n_obs = sum(!is.na(y)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  none <- out$item[out$n_left + out$n_right == 0]
  if (length(none)) message("censor: no values at or beyond the limits of ", paste(none, collapse = ", "))
  out
}

# the censoring part of S (rtmb_build): indicator matrices of left / right censored values, the
# values replaced by the limit (the density term of a censored value is multiplied by 0)
censor_build <- function(S, cs) {
  if (is.null(cs)) return(S)
  N <- S$N; nC <- length(S$C)
  S$cens <- data.frame(j = match(cs$item, S$C), lower = cs$lower, upper = cs$upper)
  S$Lc <- S$Rc <- matrix(0, N, nC)
  for (r in seq_len(nrow(S$cens))) {
    j <- S$cens$j[r]; ok <- S$Mc[, j] == 1; y <- S$Yc[, j]
    S$Lc[, j] <- 1 * (ok & y <= S$cens$lower[r]); S$Rc[, j] <- 1 * (ok & y >= S$cens$upper[r])
    S$Yc[S$Lc[, j] == 1, j] <- S$cens$lower[r]; S$Yc[S$Rc[, j] == 1, j] <- S$cens$upper[r]
  }
  S
}

# log F(lower) and log S(upper) of continuous indicator j at the location mu, log scale lsc
# (and SHASH skewness e, tail weight dl; e = NULL for normal residuals)
cens_logp <- function(lim, mu, lsc, e, dl, upper) {
  x <- (lim - mu) / exp(lsc)
  u <- if (is.null(e)) x else T_inv(x, e, dl)
  RTMB::pnorm(u, lower.tail = !upper, log.p = TRUE)               # RTMB generic: also for AD types
}

# the header line of a censored fit
censor_txt <- function(cs) {
  if (is.null(cs)) return(NULL)
  lim <- function(v, raw, lg) if (!is.finite(v)) "" else if (lg) sprintf("%s (log %s)", format(signif(v, 4)), format(raw)) else format(signif(v, 4))
  it <- vapply(seq_len(nrow(cs)), function(i) {
    r <- cs[i, ]
    paste0(r$item, ": ", paste(c(if (is.finite(r$lower)) sprintf("<= %s %d left", lim(r$lower, r$lower_raw, r$log), r$n_left),
                                 if (is.finite(r$upper)) sprintf(">= %s %d right", lim(r$upper, r$upper_raw, r$log), r$n_right)), collapse = ", "))
  }, "")
  sprintf("  Censored (only beyond the limit is known): %s%s", paste(utils::head(it, 3), collapse = "; "),
          if (length(it) > 3) sprintf("; ... (%d indicators, %d censored values)", length(it), sum(cs$n_left + cs$n_right)) else "")
}
