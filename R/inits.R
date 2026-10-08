# Default initial values: measurement loadings 1, cross-loadings 0 (spike-and-slab:
# included, slab 0), intercepts and residual scales from the data, factor SDs 0.5,
# ALD mixing variables at their scale; small jitter per chain. Starting
# near the simple-structure solution matters: with free cross-loadings on every
# indicator and weak priors, a factor can rotate into a separate posterior mode.

make_inits <- function(spec, data, code, seed) {
  force(spec); force(data); force(code); force(seed)
  pt <- code$partab; N <- nrow(data)
  sv <- start_values(spec, data)
  function(ch) {
    set.seed(seed + 100 * ch)                # runs in a JAGS worker process
    jit <- function(v, rel = FALSE) if (rel) v * (1 + stats::runif(1, -0.1, 0.1)) else v + stats::rnorm(1, 0, 0.05)
    ini <- list()
    # the scale of each latent variable (loadings x c, persons / c) starts at a different place in
    # every chain, so that R-hat can see a slowly mixing scale
    cs <- stats::setNames(exp(stats::runif(length(spec$factors), -0.3, 0.3)), spec$factors)
    for (f in names(sv$person)) ini[[f]] <- sv$person[[f]] / cs[[f]]
    q_init <- function(sx, s0) ini[[sprintf("e_%s", sx)]] <<- rep(s0, N)   # ALD mixing variables
    for (i in seq_len(nrow(pt))) {
      r <- pt[i, ]
      if (is.na(r$node)) next
      if (r$kind == "loading") {
        first <- spec$ld$first[spec$ld$lhs == r$lhs & spec$ld$rhs == r$rhs][1]
        v0 <- if (!isTRUE(first)) 0 else if (!is.null(l0 <- sv$load[[r$lhs]][r$rhs]) && !is.na(l0)) unname(l0) * cs[[r$lhs]] else 1
        if (!is.na(r$incl)) { ini[[r$incl]] <- 1; ini[[sprintf("slab_%s", r$node)]] <- jit(v0) }
        else if (is.null(ini[[r$node]])) ini[[r$node]] <- jit(v0, rel = isTRUE(first))     # positive loadings: multiplicative jitter
      } else if (r$kind == "beta") {
        if (!is.na(r$incl)) { ini[[r$incl]] <- 1; ini[[sprintf("slab_%s", r$node)]] <- jit(0) }
        else ini[[r$node]] <- jit(0)
      } else if (r$kind == "intercept") {
        v <- data[[r$lhs]]
        ini[[r$node]] <- if (r$lhs %in% spec$binary) jit(stats::qlogis(min(max(mean(v, na.rm = TRUE), 0.02), 0.98)))
                         else jit(mean(v, na.rm = TRUE))
      } else if (r$kind %in% c("resid", "resid_ald")) {
        s0 <- 0.8 * stats::sd(data[[r$lhs]], na.rm = TRUE)
        f <- spec$family$ind[[r$lhs]]
        if (is_q(f)) { ini[[sprintf("sdr_%s", safe(r$lhs))]] <- jit(s0, TRUE); q_init(safe(r$lhs), s0 / sqrt(ald_var(f$p))) }
        else ini[[r$node]] <- jit(s0, TRUE)
      } else if (r$kind == "fsd") ini[[r$node]] <- jit(sv$sd[[r$lhs]] %||% 0.5, TRUE)
      else if (r$kind == "corr") ini[[r$node]] <- jit(0)
      else if (r$kind %in% c("alpha", "delta", "psi", "kappa")) {          # moderation effects start near 0
        if (!is.na(r$incl)) { ini[[r$incl]] <- 1; ini[[sprintf("slab_%s", r$node)]] <- stats::rnorm(1, 0, 0.02) }
        else if (is.null(ini[[r$node]])) ini[[r$node]] <- stats::rnorm(1, 0, 0.02)
      }
    }
    for (f in spec$factors) {
      ff <- spec$family$factor[[f]]
      if (!is_q(ff)) next
      s0 <- (if (is.na(spec$fsd[[f]])) 0.5 else spec$fsd[[f]]) / sqrt(ald_var(ff$p))
      if (is.na(spec$fsd[[f]])) ini[[sprintf("sd_%s", f)]] <- jit(0.5, TRUE)
      q_init(f, s0)
    }
    for (f in spec$factors) if (grepl(sprintf("mu_load_%s ~", f), code$code, fixed = TRUE)) {
      l0 <- sv$load[[f]]
      ini[[sprintf("mu_load_%s", f)]] <- if (length(l0) && any(!is.na(l0))) mean(log(l0), na.rm = TRUE) + log(cs[[f]]) else 0
    }
    for (v in spec$ordinal %||% character()) {                   # thresholds from the cumulative proportions
      y <- spec$ord$Y[, v]; K <- spec$ord$K[match(v, spec$ord$items)]
      cp <- cumsum(tabulate(y, K))[-K] / sum(!is.na(y))
      bb <- stats::qlogis(pmin(pmax(cp, 0.01), 0.99)) / 1.2
      if (length(bb) > 1) for (k in 2:length(bb)) bb[k] <- max(bb[k], bb[k - 1] + 0.1)
      if (spec$ord$type[[v]] != "grm") {                       # partial credit: steps at the thresholds, step ratios 1
        ini[[sprintf("b_%s", safe(v))]] <- bb + stats::rnorm(K - 1, 0, 0.05)
        if (spec$ord$type[[v]] == "tppcm" && K > 2) ini[[sprintf("r_%s", safe(v))]] <- c(NA, rep(1, K - 2))
        next
      }
      ini[[sprintf("b_%s", safe(v))]] <- c(bb[1] + stats::rnorm(1, 0, 0.05), rep(NA, K - 2))
      if (K > 2) ini[[sprintf("inc_%s", safe(v))]] <- c(NA, diff(bb))
    }
    for (g in unique(regmatches(code$code, gregexpr("(?<=sigma_slab_)[A-Za-z0-9_]+(?= ~)", code$code, perl = TRUE))[[1]])) {
      if (grepl(sprintf("p_incl_%s ~", g), code$code, fixed = TRUE)) ini[[sprintf("p_incl_%s", g)]] <- 0.5
      ini[[sprintf("sigma_slab_%s", g)]] <- 0.2
    }
    ini
  }
}

# Starting values from the data. For each latent variable, a proxy score is the mean of its
# standardized first-line indicators (signed by fixed loadings); its reliability (alpha) corrects
# the attenuation. Free loadings of a unit-variance latent variable start at the implied slopes
# (continuous: cov(x, f); binary: logistic slope), persons at the shrunken proxy, and a free SD
# (all loadings fixed) at the SD of the proxy in the units of the indicators. Starting the loadings
# at 1 instead puts a latent variable measured by many indicators on the wrong scale, and the
# Gibbs sampler moves along the scale (loadings x c, persons / c) very slowly.
start_values <- function(spec, data) {
  ld <- spec$ld[!spec$ld$zero, ]; out <- list(load = list(), person = list(), sd = list())
  for (f in spec$factors) {
    r <- ld[ld$lhs == f & ld$first, ]
    if (nrow(r) < 2) next
    fx <- suppressWarnings(as.numeric(r$fixed)); sg <- ifelse(is.na(fx), 1, sign(fx))
    X <- as.matrix(data[r$rhs]); storage.mode(X) <- "double"
    mu <- colMeans(X, na.rm = TRUE); sdx <- apply(X, 2, stats::sd, na.rm = TRUE)
    Z <- sweep(sweep(X, 2, mu), 2, sdx / sg, `/`)
    C <- stats::cor(Z, use = "pairwise.complete.obs"); k <- ncol(Z); mc <- mean(C[upper.tri(C)], na.rm = TRUE)
    rel <- min(max(k * mc / (1 + (k - 1) * mc), 0.3), 0.95)
    s <- rowMeans(Z, na.rm = TRUE); s[!is.finite(s)] <- 0; s <- (s - mean(s)) / stats::sd(s)
    if (is.na(spec$fsd[[f]])) {                                           # loadings fixed, SD free
      if (all(!is.na(fx))) {
        u <- rowMeans(sweep(sweep(X, 2, mu), 2, fx, `/`), na.rm = TRUE); u[!is.finite(u)] <- 0
        out$sd[[f]] <- max(stats::sd(u) * sqrt(rel), 1e-3)
        out$person[[f]] <- (u - mean(u)) * rel
      }
      next
    }
    alpha <- function(k) min(max(k * mc / (1 + (k - 1) * mc), 0.3), 0.95)
    l0 <- vapply(seq_len(k), function(i) {                              # rest score: without the item itself
      x <- X[, i]; p <- rowMeans(Z[, -i, drop = FALSE], na.rm = TRUE); ok <- !is.na(x) & is.finite(p)
      p <- (p - mean(p[ok])) / stats::sd(p[ok])
      b <- if (r$rhs[i] %in% spec$binary) {
        g <- tryCatch(suppressWarnings(stats::glm.fit(cbind(1, p[ok]), x[ok], family = stats::binomial())), error = function(e) NULL)
        if (is.null(g)) NA_real_ else g$coefficients[2]
      } else stats::cov(x[ok], p[ok])
      b / sqrt(alpha(k - 1)) / spec$fsd[[f]]
    }, 0)
    l0 <- ifelse(is.na(fx), pmin(pmax(l0, 0.05), 4), NA_real_)
    out$load[[f]] <- stats::setNames(l0, r$rhs)
    out$person[[f]] <- s * sqrt(rel) * spec$fsd[[f]]
  }
  out
}
