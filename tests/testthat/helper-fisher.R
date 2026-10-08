# Numerical expected Fisher information about (z1, z2) at one value theta of the first latent
# variable (speed at z2 = 0, covariates / moderators at xb), from the log-likelihood code itself
# (ind_ll()): minus the expected Hessian, by central differences, of each indicator's log-likelihood,
# the expectation over its outcomes (binary, ordinal categories, Gauss-Hermite or Gauss-Legendre
# nodes of a continuous residual, plus the censored masses). `limits`: optional named list of
# c(lower, upper) per continuous indicator (data scale) that replaces the censoring of the fit.
# Returns the 2 x 2 information summed over each group of indicators (cat, rt), as cr_zero() lists.
fisher_num <- function(fit, theta, xb = NULL, limits = NULL, h = 2e-4, k = 60) {
  S <- fit$S; p <- fit$rtmb$par; nF <- S$nF
  xb <- xb %||% if (ncol(S$X)) colMeans(S$X) else numeric(0)
  st <- birt:::cr_setup(p, S, theta, xb)
  f1i <- S$f[[1]]; mu1 <- birt:::lin_x(p$beta1, f1i$cols, matrix(xb, 1, dimnames = list(NULL, colnames(S$X))))[1]
  z1 <- (theta - mu1) / st$s1
  one <- function(type, j, y, w, cens = NULL) {           # outcomes y (M) with probabilities w
    M <- length(y); S1 <- S; S1$N <- M
    S1$X <- matrix(xb, M, length(xb), byrow = TRUE, dimnames = list(NULL, colnames(S$X)))
    S1$Mb <- matrix(0, M, length(S$B)); S1$Yb <- S1$Mb; S1$Mc <- matrix(0, M, length(S$C)); S1$Yc <- S1$Mc
    S1$Lc <- S1$Rc <- S1$Mc
    if (length(S$O)) { S1$Mo <- matrix(0, M, length(S$O)); S1$ind <- lapply(seq_along(S$O), function(i) matrix(0, M, S$K[i])) }
    if (type == "b") { S1$Mb[, j] <- 1; S1$Yb[, j] <- y }
    if (type == "o") { S1$Mo[, j] <- 1; S1$ind[[j]][cbind(seq_len(M), y)] <- 1 }
    if (type == "c") {
      S1$Mc[, j] <- 1; S1$Yc[, j] <- y
      if (!is.null(cens)) { S1$cens <- data.frame(j = j, lower = cens[1], upper = cens[2]); S1$Lc[, j] <- 1 * (y <= cens[1]); S1$Rc[, j] <- 1 * (y >= cens[2]) }
      else S1$cens <- NULL
    }
    rows <- seq_len(M)
    ll <- function(a, b) birt:::ind_ll(p, S1, birt:::lat_values(p, S1, rep(a, M), rep(b, M), rows), rows)
    a <- z1; b <- 0
    if (is.null(w)) w <- exp(ll(a, b))                     # categorical: the probabilities of the outcomes
    h11 <- (ll(a + h, b) - 2 * ll(a, b) + ll(a - h, b)) / h^2
    if (nF == 1) return(list(i11 = -sum(w * h11), i12 = 0, i22 = 0))
    h22 <- (ll(a, b + h) - 2 * ll(a, b) + ll(a, b - h)) / h^2
    h12 <- (ll(a + h, b + h) - ll(a + h, b - h) - ll(a - h, b + h) + ll(a - h, b - h)) / (4 * h^2)
    list(i11 = -sum(w * h11), i12 = -sum(w * h12), i22 = -sum(w * h22))
  }
  plus <- function(A, B) Map(`+`, A, B)
  cat <- list(i11 = 0, i12 = 0, i22 = 0); rt <- cat
  for (j in seq_along(S$B)) cat <- plus(cat, one("b", j, c(0, 1), NULL))
  for (j in seq_along(S$O)) cat <- plus(cat, one("o", j, seq_len(S$K[j]), NULL))
  gh <- statmod::gauss.quad.prob(k, "normal"); gl <- statmod::gauss.quad(k, "legendre")
  for (j in seq_along(S$C)) {
    mu0 <- st$P0$mu[1, j]; sd0 <- exp(st$P0$lsc[1, j]); s <- match(j, S$sh); e0 <- st$P0$e[1, j]; dl <- if (is.na(s)) 1 else exp(p$ldl[s])
    Tz <- function(z) if (is.na(s)) z else birt:::T_std(z, e0, dl)
    Ti <- function(x) if (is.na(s)) x else birt:::T_inv(x, e0, dl)
    lim <- limits[[S$C[j]]] %||% (if (!is.null(S$cens) && j %in% S$cens$j) unlist(S$cens[S$cens$j == j, c("lower", "upper")]))
    if (is.null(lim)) { rt <- plus(rt, one("c", j, mu0 + sd0 * Tz(gh$nodes), gh$weights)); next }
    zl <- if (is.finite(lim[1])) Ti((lim[1] - mu0) / sd0) else -10; zu <- if (is.finite(lim[2])) Ti((lim[2] - mu0) / sd0) else 10
    zq <- (zu - zl) / 2 * gl$nodes + (zu + zl) / 2; wq <- (zu - zl) / 2 * gl$weights * stats::dnorm(zq)
    yy <- c(mu0 + sd0 * Tz(zq), if (is.finite(lim[1])) lim[1], if (is.finite(lim[2])) lim[2])
    ww <- c(wq, if (is.finite(lim[1])) stats::pnorm(zl), if (is.finite(lim[2])) stats::pnorm(zu, lower.tail = FALSE))
    rt <- plus(rt, one("c", j, yy, ww, cens = unname(lim)))
  }
  list(cat = cat, rt = rt, all = plus(cat, rt), nF = nF, s1 = st$s1,
       info = function(I) birt:::cr_rel(I, nF, st$s1)$info)
}
