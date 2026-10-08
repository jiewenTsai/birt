#' Quantile fits: several quantiles of one latent variable or indicator
#'
#' `birt(..., family = list(speed = ald(c(.1, .5, .9))))` fits one model per quantile
#' and returns a `birt_quantile` object: `fits` (list of `birt` objects), `p`, `target`.
#' @name birt_quantile
#' @param object,x A `birt_quantile` object.
#' @param ... Unused.
NULL

#' @rdname birt_quantile
#' @export
estimates.birt_quantile <- function(object, ...)
  do.call(rbind, Map(function(f, p) cbind(p = p, estimates(f)), object$fits, object$p))

#' @rdname birt_quantile
#' @export
coef.birt_quantile <- function(object, ...) {
  cf <- lapply(object$fits, stats::coef)
  nm <- unique(unlist(lapply(cf, names)))
  out <- vapply(cf, function(v) v[nm], numeric(length(nm)))
  dimnames(out) <- list(nm, sprintf("p=%s", format(object$p)))
  out
}

#' @rdname birt_quantile
#' @export
reliability.birt_quantile <- function(object, ...) {
  out <- t(vapply(object$fits, function(f) f$reliability, object$fits[[1]]$reliability))
  rownames(out) <- sprintf("p=%s", format(object$p)); out
}

#' @rdname birt_quantile
#' @export
fit_indices.birt_quantile <- function(object, ...)
  do.call(rbind, Map(function(f, p) cbind(p = p, f$fit_indices), object$fits, object$p))

#' @rdname birt_quantile
#' @export
print.birt_quantile <- function(x, digits = 3, ...) {
  cat(sprintf("birt quantile fits for '%s' at p = %s\n", x$target, paste(format(x$p), collapse = ", ")))
  e <- estimates(x)
  rg <- e[e$op == "~" & (e$lhs == x$target | x$target == ".continuous"), ]
  if (nrow(rg)) {
    use_adj <- any(!is.na(rg$sd_adj))
    if (use_adj) { rg$q025 <- ifelse(is.na(rg$q025_adj), rg$q025, rg$q025_adj); rg$q975 <- ifelse(is.na(rg$q975_adj), rg$q975, rg$q975_adj) }
    cat(sprintf("Regression coefficients: mean [95%% %s]\n", if (use_adj) "interval, adjusted for the ALD working likelihood (Yang, Wang & He, 2016)" else "interval"))
    key <- unique(paste(rg$lhs, "~", rg$rhs))
    tab <- data.frame(term = key)
    for (p in x$p) {
      r <- rg[rg$p == p, ]; i <- match(key, paste(r$lhs, "~", r$rhs))
      tab[[sprintf("p=%s", format(p))]] <- sprintf("%.*f [%.*f, %.*f]", digits, r$est[i], digits, r$q025[i], digits, r$q975[i])
    }
    print_table(tab, digits)
  }
  cat("Reliability\n"); print_table(data.frame(p = x$p, reliability(x), check.names = FALSE), digits)
  fi <- tryCatch(fit_indices(x), error = function(e) NULL)
  if (!is.null(fi) && nrow(fi) && all(c("DIC", "LOOIC") %in% names(fi))) {
    cat("Fit indices\n"); print_table(fi[, c("p", "DIC", "WAIC", "LOOIC", "SE_LOOIC")], digits)
  }
  invisible(x)
}

#' @rdname birt_quantile
#' @export
summary.birt_quantile <- function(object, ...) { print(object, ...); invisible(object) }

#' @rdname birt_quantile
#' @param digits Number of decimals.
#' @param terms Regression terms to plot (default: all regressions of the target).
#' @export
plot.birt_quantile <- function(x, terms = NULL, ...) {
  e <- estimates(x); rg <- e[e$op == "~" & e$lhs == x$target, ]
  key <- unique(paste(rg$lhs, "~", rg$rhs)); if (!is.null(terms)) key <- intersect(key, terms)
  if (!length(key)) { message("no regression coefficients of ", x$target, " to plot"); return(invisible(NULL)) }
  nc <- ceiling(sqrt(length(key)))
  op <- graphics::par(mfrow = c(ceiling(length(key) / nc), nc)); on.exit(graphics::par(op))
  for (k in key) {
    r <- rg[paste(rg$lhs, "~", rg$rhs) == k, ]; r <- r[order(r$p), ]
    graphics::plot(r$p, r$est, type = "n", ylim = range(r$q025, r$q975, 0), xlab = "quantile p", ylab = "coefficient", main = k, ...)
    graphics::polygon(c(r$p, rev(r$p)), c(r$q025, rev(r$q975)), col = "grey85", border = NA)
    graphics::lines(r$p, r$est, type = "b", pch = 19); graphics::abline(h = 0, lty = 3)
  }
  invisible(rg)
}

# ---- adjusted posterior intervals under the ALD working likelihood ----------------------
# Yang, Wang & He (2016): with an ALD(tau, sigma) working likelihood the posterior covariance
# of quantile regression coefficients is sigma / n * D1^-1 rather than their sampling
# covariance tau (1 - tau) / n * D1^-1 D0 D1^-1 (D0 = E[x x'], D1 = E[f(0 | x) x x']). With
# latent variables the posterior covariance Sigma also carries measurement uncertainty, so
# only the likelihood part is replaced:
#   V = Sigma - sigma / (n f) D0^-1 + tau (1 - tau) / (n f^2) D0^-1,
# with f the density of the residuals at 0 (kernel estimate, iid form as in quantreg) and D0
# averaged over posterior draws of the person values. For an observed outcome V equals the
# adjusted covariance of Yang et al. Applied to the linear coefficients of every ALD target
# (the latent regression of an ALD latent variable; the intercept and loadings of an ALD
# indicator).
ald_adjust <- function(fit) {
  e <- fit$estimates; sp <- fit$spec
  e$sd_adj <- e$q025_adj <- e$q975_adj <- NA_real_
  tg <- c(names(Filter(is_q, sp$family$factor)), names(Filter(is_q, sp$family$ind)))
  if (!length(tg)) return(e)
  M <- as.matrix(fit$draws); N <- fit$data$N
  sub <- unique(round(seq(1, nrow(M), length.out = min(nrow(M), 200))))
  value <- function(x, d) {                                  # regressor x in draw d
    if (x == "(Intercept)") rep(1, N)
    else if (x %in% sp$factors) M[d, sprintf("%s[%d]", x, seq_len(N))]
    else fit$data[[x]]
  }
  outcome <- function(t, d) if (t %in% sp$factors) M[d, sprintf("%s[%d]", t, seq_len(N))] else fit$data[[t]]
  for (t in tg) {
    lat <- t %in% sp$factors
    fam <- if (lat) sp$family$factor[[t]] else sp$family$ind[[t]]
    if (lat) {
      rows <- which(e$op == "~" & e$lhs == t & e$kind == "beta" & !e$fixed)
      regs <- e$rhs[rows]
      allr <- which(e$op == "~" & e$lhs == t & e$kind == "beta")      # incl. fixed terms
      sn <- sprintf("scale_%s", t)
    } else {
      i0 <- which(e$op == "~1" & e$lhs == t & !e$fixed); il <- which(e$op == "=~" & e$rhs == t & !e$fixed)
      rows <- c(i0, il); regs <- c(rep("(Intercept)", length(i0)), e$lhs[il])
      allr <- c(which(e$op == "~1" & e$lhs == t), which(e$op == "=~" & e$rhs == t))
      sn <- sprintf("sigma_%s", safe(t))
    }
    if (!length(rows) || any(is.na(e$node[rows]))) next
    allregs <- if (lat) e$rhs[allr] else c("(Intercept)", e$lhs[allr[-1]])
    sig <- if (sn %in% colnames(M)) mean(M[, sn]) else {
      sdf <- sp$fsd[[t]]; if (is.null(sdf) || is.na(sdf)) next; sdf / sqrt(ald_var(fam$p))
    }
    D0 <- 0; fh <- numeric(length(sub)); k <- 0
    for (d in sub) {
      k <- k + 1
      X <- vapply(regs, value, numeric(N), d = d); D0 <- D0 + crossprod(X) / N
      # residual at the p-quantile: outcome - all linear terms (fixed loadings included)
      cf <- vapply(allr, function(r) if (is.na(e$node[r])) e$est[r] else M[d, e$node[r]], 0)
      Xa <- vapply(allregs, value, numeric(N), d = d)
      u <- outcome(t, d) - as.vector(Xa %*% cf)
      dd <- stats::density(u, bw = "nrd0", n = 512, from = -3 * stats::sd(u), to = 3 * stats::sd(u))
      fh[k] <- stats::approx(dd$x, dd$y, 0)$y
    }
    D0 <- D0 / length(sub); f0 <- mean(fh)
    S <- stats::cov(M[, e$node[rows], drop = FALSE]); Di <- solve(D0)
    V <- S + (fam$p * (1 - fam$p) / f0^2 - sig / f0) / N * Di
    sd <- sqrt(pmax(diag(V), 0))
    e$sd_adj[rows] <- sd; e$q025_adj[rows] <- e$est[rows] - 1.959964 * sd; e$q975_adj[rows] <- e$est[rows] + 1.959964 * sd
  }
  e
}
