# Marginal log-likelihood per person and draw ---------------------------------------
#
# The first latent variable is integrated by Gauss-Hermite quadrature (normal) or by
# Gauss-Laguerre over the ALD mixing variable e and Gauss-Hermite given e (ALD). Given
# it, the second latent variable is integrated analytically when all continuous
# indicators are normal (their covariance is v l l' + D; Woodbury identity and the matrix
# determinant lemma), with an extra Gauss-Laguerre layer when the second variable is ALD;
# with ALD indicators it is integrated by quadrature as well. Binary indicators must load
# on the first latent variable only. Missing values are skipped.

# Nodes (N x G) and log weights (G) for a latent variable with location m (N-vector).
# Normal: Gauss-Hermite; ALD: Gauss-Laguerre over e x Gauss-Hermite given e; an integrand
# with kinks (ALD indicators): a dense grid on the closed-form
# density (trapezoid rule).
latent_nodes <- function(m, Fi, s, Q, kinked = FALSE, G = 121) {
  fam <- Fi$fam$family; sc <- Fi$sc[s]
  if (fam == "normal" && !kinked) {
    gh <- statmod::gauss.quad.prob(Q, "normal")
    return(list(X = outer(m, sc * gh$nodes, "+"), lw = log(gh$weights)))
  }
  if (fam == "ald" && !kinked) {
    gl <- statmod::gauss.quad(10, "laguerre"); gh <- statmod::gauss.quad.prob(15, "normal")
    k <- ald_k(Fi$fam$p); e <- sc * gl$nodes
    shift <- rep(k[1] * e, each = 15); sdv <- rep(sqrt(k[2] * sc * e), each = 15)
    return(list(X = outer(m, shift + sdv * rep(gh$nodes, 10), "+"),
                lw = rep(log(gl$weights), each = 15) + rep(log(gh$weights), 10)))
  }
  if (fam == "normal") { ctr <- 0; sdv <- sc; ld <- function(x) stats::dnorm(x, 0, sc, log = TRUE) }
  else { p <- Fi$fam$p; ctr <- sc * ald_k(p)[1]; sdv <- sc * sqrt(ald_var(p)); ld <- function(x) dald_log(x, 0, sc, p) }
  grid <- ctr + sdv * seq(-7, 7, length.out = G)
  lw <- ld(grid); lw <- lw - (max(lw) + log(sum(exp(lw - max(lw)))))
  list(X = outer(m, grid, "+"), lw = lw)
}

# posterior mean and SD of each person's sampled f1 (NULL when f1 was not monitored)
f1_moments <- function(P) {
  th_cols <- sprintf("%s[%d]", P$f1, seq_len(P$N))
  if (!all(th_cols %in% colnames(P$M))) return(list(m = NULL, s = NULL))
  D <- P$M[, th_cols, drop = FALSE]; list(m = colMeans(D), s = pmax(apply(D, 2, stats::sd), 1e-3))
}

loglik_compute <- function(P, obs, Q = 41, QL = 10, G1 = 161, G2 = 121) {
  if (any(P$a2 != 0) || any(P$o2 != 0)) stop("the marginal likelihood needs categorical indicators to load on the first latent variable only")
  N <- P$N
  gl <- statmod::gauss.quad(QL, "laguerre")
  Y <- obs$Y; Tm <- obs$T; Yo <- P$Yo; O <- if (is.null(Yo)) character() else colnames(Yo)
  oY <- !is.na(Y); Y0 <- ifelse(oY, Y, 0); oT <- !is.na(Tm); T0 <- ifelse(oT, Tm, 0)
  C <- ncol(Tm); nobs <- rowSums(oT); B <- colnames(Y); Cn <- colnames(Tm)
  is_ald <- !is.na(P$pc); any_ald <- any(is_ald)
  softplus <- function(x) ifelse(x > 30, x, log1p(exp(x)))
  has2 <- !is.null(P$F2) && C > 0 && any(P$l2 != 0)
  # kinks in the integrand call for grids instead of Gauss rules
  kink1 <- any(is_ald & colSums(abs(P$l1)) > 0)
  kink2 <- any(is_ald & colSums(abs(P$l2)) > 0)
  LL <- matrix(NA_real_, P$S, N)
  mom <- P$f1_mom %||% f1_moments(P); pm0 <- mom$m; ps0 <- mom$s          # sampled f1: centres of the adaptive rule
  for (s in seq_len(P$S)) {
    # loadings, shifts and residual SDs per person (moderation)
    cols <- function(vars, f) matrix(vapply(vars, function(v) rep_len(f(v), N), numeric(N)), N, length(vars))
    A1 <- cols(B, function(v) P$a1[s, v] * mod_mult(P, s, paste(P$f1, v)))
    DB <- cols(B, function(v) mod_shift(P, s, v))
    L1 <- cols(Cn, function(v) P$l1[s, v] * mod_mult(P, s, paste(P$f1, v)))
    L2 <- cols(Cn, function(v) if (is.na(P$f2)) 0 else P$l2[s, v] * mod_mult(P, s, paste(P$f2, v)))
    DC <- cols(Cn, function(v) mod_shift(P, s, v))
    SG <- cols(Cn, function(v) P$sig[s, match(v, Cn)] * mod_sd(P, s, v))
    dens <- function(mu) {                         # sum of log densities of the RTs, N-vector
      out <- rep(0, N)
      for (c in seq_len(C)) {
        ld <- switch(P$ct[c], normal = stats::dnorm(T0[, c], mu[, c], SG[, c], log = TRUE),
                     ald = dald_log(T0[, c], mu[, c], SG[, c], P$pc[c]))
        out <- out + ifelse(oT[, c], ld, 0)
      }
      out
    }
    if (has2 && !any_ald) {
      Pm <- oT / SG^2
      A <- rowSums(Pm * L2^2)
      ldD <- rowSums(oT * log(SG^2))
      analytic <- function(mu0, m2, v) {             # N(m2, v) for the second factor
        den <- 1 + v * A
        R <- (T0 - mu0 - L2 * m2) * oT
        q2 <- rowSums(R * Pm * L2)
        -0.5 * nobs * log(2 * pi) - 0.5 * (ldD + log(den)) - 0.5 * (rowSums(R^2 * Pm) - v * q2^2 / den)
      }
    }
    sd1 <- mod_sd(P, s, P$f1)
    m1 <- mean_cov(P$F1, s, N)
    n1 <- latent_nodes(m1, P$F1, s, Q, kink1, G1)
    if (!identical(sd1, 1)) n1$X <- m1 + (n1$X - m1) * sd1         # V(f1) ~ z: the nodes scale with the person's SD
    if (has2) m2c <- mean_cov(P$F2, s, N)
    # ordinal indicators: loadings, shifts and thresholds at draw s
    oo <- lapply(O, function(v) list(a = P$o1[s, v] * mod_mult(P, s, paste(P$f1, v)), d = mod_shift(P, s, v), b = P$thr[[v]][s, ],
                                     r = P$rat[[v]][s, ], type = P$otype[[v]], y = Yo[, v], ok = !is.na(Yo[, v])))
    lik_at <- function(th) {                                          # log density of the indicators given f1 = th
      ly <- 0
      if (ncol(Y)) {
        eta <- matrix(P$d[s, ], N, ncol(Y), byrow = TRUE) + A1 * (th - DB)
        ly <- rowSums(oY * (Y0 * eta - softplus(eta)))
      }
      for (o in oo) {                                                 # graded response or partial credit (ordinal.R)
        pr <- cat_probs(o$type, o$a, th - o$d, o$b, o$r)[cbind(seq_len(N), ifelse(o$ok, o$y, 1))]
        ly <- ly + ifelse(o$ok, log(pmax(pr, 1e-300)), 0)
      }
      lt <- 0
      if (C) {
        mu0 <- matrix(P$xi[s, ], N, C, byrow = TRUE) + L1 * th + DC
        if (!has2) lt <- dens(mu0) else {
          m2 <- m2c + mean_path(P$F2, s, th)
          F2 <- P$F2; sd2 <- F2$sc[s] * mod_sd(P, s, P$f2); s1 <- P$F1$sc[s] * sd1
          if (F2$fam$family == "normal") {                          # conditional normal (with correlation)
            if (P$r[s] != 0) m2 <- m2 + P$r[s] * sd2 * (th - m1) / s1
            v <- sd2^2 * (1 - P$r[s]^2)
            if (!any_ald) lt <- analytic(mu0, m2, v) else {
              n2 <- latent_nodes(m2, list(fam = normal(), sc = rep(1, s)), s, Q, kink2, G2)
              n2$X <- m2 + (n2$X - m2) * sqrt(v)
              lt <- lse(vapply(seq_len(ncol(n2$X)), function(q)
                n2$lw[q] + dens(mu0 + L2 * n2$X[, q]), numeric(N)))
            }
          } else {                                                  # ALD second factor
            k <- ald_k(F2$fam$p); e <- F2$sc[s] * gl$nodes
            lt <- if (!any_ald) lse(vapply(seq_along(e), function(i)
                    log(gl$weights[i]) + analytic(mu0, m2 + k[1] * e[i], k[2] * F2$sc[s] * e[i]), numeric(N)))
                  else {
                    n2 <- latent_nodes(m2, F2, s, Q, kink2, G2)
                    lse(vapply(seq_len(ncol(n2$X)), function(q)
                      n2$lw[q] + dens(mu0 + L2 * n2$X[, q]), numeric(N)))
                  }
          }
        }
      }
      ly + lt
    }
    if (P$F1$fam$family == "normal" && !kink1) {
      # adaptive Gauss-Hermite (Liu & Pierce, 1994): nodes c + s z_q with weights w_q phi(x; m1, s1) s / phi(z_q),
      # first centred at the posterior mean and SD of the person's sampled f1, then recentred at the
      # posterior moments of the previous round until the log-likelihood is stable (at most 5 rounds)
      gh <- statmod::gauss.quad.prob(Q, "normal"); s1 <- P$F1$sc[s] * sd1
      cc <- rep_len(if (is.null(pm0)) m1 else pm0, N); ss <- rep_len(if (is.null(ps0)) s1 else ps0, N); prev <- Inf
      for (round in 1:5) {
        X <- outer(cc, rep(1, Q)) + outer(ss, gh$nodes)
        LW <- vapply(seq_len(Q), function(q) log(gh$weights[q]) + stats::dnorm(X[, q], m1, s1, log = TRUE) + log(ss) -
                       stats::dnorm(gh$nodes[q], log = TRUE) + lik_at(X[, q]), numeric(N))
        ll <- lse(LW)
        if (max(abs(ll - prev)) < 1e-7) break
        W <- exp(LW - ll); cc <- rowSums(W * X); ss <- sqrt(pmax(rowSums(W * X^2) - cc^2, 1e-10)); prev <- ll
      }
      LL[s, ] <- ll
      next
    }
    LW <- vapply(seq_len(ncol(n1$X)), function(g) n1$lw[g] + lik_at(n1$X[, g]), numeric(N))
    LL[s, ] <- lse(LW)
  }
  LL
}

obs_matrices <- function(fit) {
  sp <- fit$spec; N <- fit$data$N
  mk <- function(v) if (length(v)) matrix(unlist(fit$data[v]), N, length(v), dimnames = list(NULL, v)) else matrix(0, N, 0)
  list(Y = mk(sp$binary), T = mk(sp$cont))
}

# pointwise log-likelihood, DIC, WAIC and PSIS-LOO from it
fit_indices_compute <- function(P, obs, chain, Q = 41) {
  LL <- loglik_compute(P, obs, Q = Q)
  dev <- -2 * rowSums(LL)
  r_eff <- loo::relative_eff(exp(LL), chain_id = chain)
  w <- suppressWarnings(loo::waic(LL))
  l <- suppressWarnings(loo::loo(LL, r_eff = r_eff))
  k <- l$diagnostics$pareto_k
  tab <- data.frame(Deviance = mean(dev), pD = stats::var(dev) / 2, DIC = mean(dev) + stats::var(dev) / 2,
                    WAIC = w$estimates["waic", "Estimate"], p_WAIC = w$estimates["p_waic", "Estimate"],
                    LOOIC = l$estimates["looic", "Estimate"], SE_LOOIC = l$estimates["looic", "SE"],
                    p_LOO = l$estimates["p_loo", "Estimate"], k_gt_0.7 = sum(k > 0.7))
  list(table = tab, loglik = LL, loo = l, waic = w)
}

# cost-based number of draws for the marginal likelihood (ALD adds quadrature layers)
n_draws_ll <- function(spec, n_draws) {
  alds <- any(vapply(spec$family$ind, is_q, logical(1)))
  kink1 <- alds
  f1 <- spec$family$factor[[1]]$family
  g1 <- if (kink1) 161 else c(normal = 41, ald = 150)[[f1]]
  g2 <- if (length(spec$factors) < 2) 1 else {
    f2 <- spec$family$factor[[2]]$family
    if (alds) 121 else c(normal = 1, ald = 10)[[f2]]
  }
  cost <- g1 * g2 / 41
  min(n_draws, max(100, round(n_draws / cost)))
}
