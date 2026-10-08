# RTMB engine -----------------------------------------------------------------------------
#
# Maximum likelihood for the same lavaan-style RT-IRT models, with the latent variables
# integrated out by adaptive Gauss-Hermite quadrature (AGHQ) or the Laplace approximation.
# Every latent variable is a transform of independent standard normals z1, z2:
#   f1 = x'beta1 + s1(x) z1,                         log s1 = lsd1 + x'alpha1
#   f2 = gam f1 + x'beta2 + s2(m) T(u),  u = r z1 + sqrt(1 - r^2) z2,  log s2 = lsd2 + m'alpha2
# with T the identity (normal) or the sinh-arcsinh transform sinh((asinh u + eps) / delta)
# (shash(); Jones & Pewsey, 2009), standardized to T(0) = 0, T'(0) = 1 (T_std), whose skewness
# eps may depend on moderators m. Continuous
# indicators have y = xi + sum lam f + sd(m) e with e normal or SHASH(eps(m), delta), log sd linear
# in m. Quantile effects are derivatives of the conditional p-quantiles, averaged over the
# model-implied distribution of the moderators.
#
# AGHQ: per person a product rule at the conditional mode of (z1, z2) with the Cholesky factor of
# the inverse Hessian (from the Laplace object). The nodes are data, so the AGHQ objective is an
# ordinary fixed-effect RTMB function; the centres are recomputed at the new estimates until the
# objective stops changing.

#' Control settings of the RTMB engine
#'
#' @param method `"elgm"` (default): approximate Bayesian inference for the model as an extended latent Gaussian
#'   model (Stringer, Brown & Stafford, 2023) with the priors of [dpriors()] and `prior()`.
#'   The latent Gaussian field W holds the person latents only; they are integrated out per
#'   person (`inner`), and the posterior of the remaining parameters is integrated by an outer
#'   adaptive quadrature: `k_hyper` nodes per dimension on the hyperparameters `hyper`, a
#'   Gaussian (Laplace) treatment of the other directions. See *Details*. `"aghq"`: maximum
#'   likelihood with adaptive Gauss-Hermite quadrature over the person latents (as GLMMadaptive
#'   or mirt; needed by [score_test()], [quantile_score()], [dif_tree()]); `"laplace"`: maximum
#'   likelihood, Laplace.
#' @param k Quadrature nodes per latent dimension (default 21 for one latent variable, 9 for two).
#' @param inner `method = "elgm"`: the integral over the person latents, `"aghq"` (default; the
#'   field is separable by person, so per-person AGHQ is cheap and more accurate) or
#'   `"laplace"` (the Laplace approximation of the original ELGM).
#' @param hyper `method = "elgm"`: parameters with quadrature in the outer layer. `"auto"`
#'   (default): the SHASH shape (and scale / skew moderators) of the second latent variable when
#'   it is `shash()`, otherwise `c("cor", "path", "sd")`. Keywords `"cor"` (latent correlation),
#'   `"path"` (latent regression path), `"sd"` (free latent SDs), `"shape"`, `"cross"`
#'   (cross-loadings), `"resid"` (residual log SDs), or names of internal parameters
#'   (`names(fit$rtmb$x)`). `character()` treats every direction as Gaussian (posterior mode +
#'   Laplace). `"auto"` uses at most 3 directions (the SHASH skewness and tail weight first, then
#'   skew and scale moderators, then correlation, path and SDs) and says which ones it treats as
#'   Gaussian: with 4 directions the 5^4 = 625 nodes took 5 times the outer time of 3 directions
#'   for posterior means within 0.034 SD and SDs within 2.1\% (a SHASH speed with skew and scale
#'   moderators, N = 1000, 20 items); name them in `hyper` for quadrature on all.
#' @param k_hyper Outer nodes per hyperparameter (default 11, 9, 5, 5, 3, 3 for 1 to 6
#'   hyperparameters). Each quadrature direction is scaled separately on both sides from the
#'   drop of the log posterior at +-2 SD (as in INLA's grid), so skewed posteriors are covered.
#' @param cond_modes `method = "elgm"`: at each outer node, move the Gaussian directions to
#'   their conditional posterior mode (default `TRUE`; Laplace over them). `FALSE` places them
#'   on the Gaussian conditional-mean line of the posterior mode, which is faster but thins the
#'   tails when the conditional modes move nonlinearly (e.g. SHASH skewness).
#' @param max_rounds,tol AGHQ: recompute the quadrature centres at most `max_rounds` times,
#'   until -2 log L changes by less than `tol`.
#' @param p Default quantile levels of [quantiles()].
#' @param trace Print the progress of the optimization.
#' @details
#' **ELGM.** With `method = "elgm"` the fit is the posterior mode of
#' log p(y | theta) + log p(theta), with log p(y | theta) from per-person AGHQ (or Laplace),
#' and the outer posterior p(theta | y) is integrated on a mixed grid: the `hyper` directions
#' come first in the Cholesky factor of the inverse Hessian, the other parameters sit at their
#' conditional posterior modes (`cond_modes`) and contribute their conditional (Laplace) variance. Posterior
#' means, SDs and 95\% intervals (from the quadrature moments, on the log scale for SDs and the
#' Fisher-z scale for correlations) are reported for every parameter and for [quantiles()];
#' EAP scores are mixtures over the outer nodes; `fit_indices()` gives the log marginal
#' likelihood (for Bayes factors in [compare()]). Item parameters are deliberately kept out of
#' W: with items in W the inner mode is a penalized joint maximum likelihood (incidental
#' parameter bias).
#' @references Stringer, A., Brown, P., & Stafford, J. (2023). Fast, scalable approximations to
#'   posterior distributions in extended latent Gaussian models. *Journal of Computational and
#'   Graphical Statistics, 32*, 84-98.
#' @return A list used by [birt()] (`engine = "rtmb"`).
#' @export
rtmb_control <- function(method = c("elgm", "aghq", "laplace"), k = NULL, max_rounds = 8, tol = 1e-3,
                         p = c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95), trace = FALSE,
                         inner = c("aghq", "laplace"), hyper = "auto", k_hyper = NULL,
                         cond_modes = TRUE) {
  if (!is.null(k_hyper) && (!is.numeric(k_hyper) || length(k_hyper) != 1 || k_hyper < 1)) stop("k_hyper must be a positive integer")
  list(method = match.arg(method), k = k, max_rounds = max_rounds, tol = tol, p = p, trace = trace,
       inner = match.arg(inner), hyper = as.character(hyper), k_hyper = k_hyper, cond_modes = isTRUE(cond_modes))
}

ld_sh <- function(x, eps, dl) {                     # SHASH(eps, delta) log density
  s <- dl * asinh(x) - eps
  s <- s / (1 + (s / 40)^8)^0.125                   # numerical guard: identity to 1e-7 for |s| < 8 (density there < exp(-7e5)), |s| < 40
  log(dl) + log(cosh(s)) - 0.5 * log(1 + x * x) - 0.5 * sinh(s)^2 - 0.9189385332046727
}
# Standardized SHASH transform T(u) = {sinh((asinh u + eps) / delta) - sinh(eps / delta)} / k,
# k = cosh(eps / delta) / delta, so that T(0) = 0 and T'(0) = 1: the location is the median and
# the scale is the scale near the median, which keeps eps and delta from trading off with the
# intercepts and the SD when the residual is close to normal.
# Numerical guard: the effective eps / delta is smoothly bounded at +-8 (a = 8 tanh(eps / (8 delta)),
# unchanged to within 1% for |eps / delta| < 1.4), since sinh(eps / delta) and cosh(eps / delta)
# cancel to garbage beyond that.
eff_a <- function(eps, dl) 8 * tanh(eps / dl / 8)
T_std <- function(u, eps, dl) { a <- eff_a(eps, dl); (sinh(asinh(u) / dl + a) - sinh(a)) * dl / cosh(a) }
ld_std <- function(x, eps, dl) {                    # log density of T_std(z), z ~ N(0, 1)
  a <- eff_a(eps, dl); k <- cosh(a) / dl
  ld_sh(x * k + sinh(a), a * dl, dl) + log(k)
}

# ---- model structure -------------------------------------------------------------------
rtmb_build <- function(sp, data) {
  if (!requireNamespace("RTMB", quietly = TRUE)) stop("engine = \"rtmb\" needs the RTMB package: install.packages(\"RTMB\")", call. = FALSE)
  if (any(sp$reg$label != "")) stop("labels on regression coefficients are not supported with engine = \"rtmb\"", call. = FALSE)
  fs <- sp$factors; nF <- length(fs); B <- sp$binary; C <- sp$cont; O <- sp$ordinal %||% character(); fam <- sp$family
  mo <- sp$mod %||% mod_empty()
  if (any(mo$type == "ssp")) stop("prior(\"ssp\") needs engine = \"jags\"; with engine = \"rtmb\" give free effects (E(y1) ~ z) and compare models by AIC / BIC", call. = FALSE)
  all_fam <- c(fam$factor, fam$ind)
  if (any(vapply(all_fam, function(f) f$family == "ald", TRUE)))
    stop("ald() is the working likelihood of one quantile (engine = \"jags\"). With engine = \"rtmb\" use shash() ",
         "(and V(y) ~ z) and get every quantile from one fit with quantiles(fit, p)", call. = FALSE)
  if (fam$factor[[fs[1]]]$family != "normal")
    stop("shash() is available for ", if (nF == 2) sprintf("the second latent variable (%s) and ", fs[2]) else "",
         "continuous indicators; ", fs[1], " is normal", call. = FALSE)
  mods <- function(f, ok, who) {
    bad <- setdiff(unique(c(f$scale, f$skew)), ok)
    if (length(bad)) stop(who, ": moderators must be among ", paste(ok, collapse = ", "), "; not: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  # covariates: regressions and covariate moderators (columns of data)
  mod_all <- unique(unlist(lapply(all_fam, function(f) c(f$scale, f$skew))))
  cov_mod <- setdiff(mod_all, fs)
  miss <- setdiff(cov_mod, names(data))
  if (length(miss)) stop("moderators not in the model or in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  Xn <- unique(c(sp$covs, cov_mod, sp$modvars))
  bad <- Xn[vapply(Xn, function(v) !is.numeric(data[[v]]) || anyNA(data[[v]]), TRUE)]
  if (length(bad)) stop("covariates must be numeric without missing values: ", paste(bad, collapse = ", "), call. = FALSE)
  X <- if (length(Xn)) as.matrix(data[Xn]) else matrix(0, nrow(data), 0)
  mods(fam$factor[[fs[1]]], Xn, fs[1])
  if (nF == 2) mods(fam$factor[[fs[2]]], c(Xn, fs[1]), fs[2])
  for (v in C) mods(fam$ind[[v]], c(Xn, fs), v)

  ld <- sp$ld[!sp$ld$zero, c("lhs", "rhs", "fixed", "label", "first", "prior")]
  ld$type <- ifelse(ld$rhs %in% B, "b", ifelse(ld$rhs %in% O, "o", "c"))
  ld$j <- ifelse(ld$type == "b", match(ld$rhs, B), ifelse(ld$type == "o", match(ld$rhs, O), match(ld$rhs, C)))
  ld$f <- match(ld$lhs, fs); ld$value <- suppressWarnings(as.numeric(ld$fixed))
  ld$kind <- ifelse(ld$fixed != "", "fixed", ifelse(ld$first, "pos", "real"))
  key <- ifelse(ld$label != "", paste0("lab:", ld$label), paste0("row:", seq_len(nrow(ld))))
  ld$kind <- ld$kind[match(key, key)]
  ld$idx <- NA_integer_
  for (k in c("pos", "real")) { i <- ld$kind == k; ld$idx[i] <- match(key[i], unique(key[i])) }

  Yb <- as.matrix(data[B]); Mb <- !is.na(Yb); Yb[!Mb] <- 0
  Yc <- as.matrix(data[C]); Mc <- !is.na(Yc); Yc[!Mc] <- 0
  storage.mode(Yb) <- storage.mode(Yc) <- "double"; storage.mode(Mb) <- storage.mode(Mc) <- "double"
  famc <- fam$ind[C]
  if (any(vapply(famc, function(f) f$family == "normal" && length(f$skew), TRUE))) stop("normal() has no skew; use shash(skew = )", call. = FALSE)
  sh <- which(vapply(famc, function(f) f$family == "shash", TRUE))
  tabmod <- function(js, what) {
    r <- do.call(rbind, lapply(js, function(j) { m <- famc[[j]][[what]]; if (length(m)) data.frame(j = j, m = m) }))
    if (is.null(r)) data.frame(j = integer(), m = character()) else r
  }
  kap <- tabmod(seq_along(C), "scale"); eta <- tabmod(sh, "skew")
  rg <- sp$reg
  f_info <- lapply(seq_len(nF), function(i) {
    f <- fs[i]; ff <- fam$factor[[f]]; r <- rg[rg$lhs == f, ]
    list(name = f, cols = match(r$x[!r$is_factor], Xn), covs = r$x[!r$is_factor],
         cov_prior = r$prior[!r$is_factor], path_prior = c(r$prior[r$is_factor], "")[1],
         path = any(r$is_factor), free_sd = is.na(sp$fsd[[f]]), sd = sp$fsd[[f]],
         shash = ff$family == "shash", scale = ff$scale %||% character(), skew = ff$skew %||% character(),
         scale_prior = ff$scale_prior %||% character())
  })
  N <- nrow(data)
  S <- list(N = N, nF = nF, fs = fs, B = B, C = C, O = O, Xn = Xn, X = X, Yb = Yb, Mb = Mb, Yc = Yc, Mc = Mc,
            ld = ld, sh = sh, kap = kap, eta = eta, f = f_info, cor = isTRUE(sp$cov_free), famc = famc)
  S$ld_b <- lapply(seq_along(B), function(j) which(ld$type == "b" & ld$j == j))
  S$ld_c <- lapply(seq_along(C), function(j) which(ld$type == "c" & ld$j == j))
  S$ld_o <- lapply(seq_along(O), function(j) which(ld$type == "o" & ld$j == j))
  if (length(O)) {                                                       # ordinal: one-hot categories, missing rows
    S$K <- sp$ord$K[match(O, sp$ord$items)]; S$otype <- unname(sp$ord$type[O])
    # parameter positions per item: grm t1 + linc (log increments), gpcm / tppcm pb (free steps), tppcm lr (log step ratios)
    o1 <- 0; o2 <- 0; o3 <- 0; o4 <- 0
    S$tix <- lapply(seq_along(O), function(j) { K <- S$K[j]; ty <- S$otype[j]
      if (ty == "grm") { o1 <<- o1 + 1; r <- list(t1 = o1, linc = o2 + seq_len(K - 2)); o2 <<- o2 + K - 2; return(r) }
      r <- list(pb = o3 + seq_len(K - 1)); o3 <<- o3 + K - 1
      if (ty == "tppcm") { r$lr <- o4 + seq_len(K - 2); o4 <<- o4 + K - 2 }
      r })
    S$ind <- lapply(seq_along(O), function(j) { y <- sp$ord$Y[, O[j]]; M <- matrix(0, N, S$K[j]); ok <- !is.na(y); M[cbind(which(ok), y[ok])] <- 1; M })
    S$Mo <- 1 - 1 * is.na(sp$ord$Y[, O, drop = FALSE])
  }
  # moderation effects: free cells (af, bf) and labelled groups (ac, bc)
  act <- mo[mo$type != "none", ]
  S$mod <- if (nrow(act)) {
    act$par <- ifelse(act$type == "common", ifelse(act$kind == "alpha", "ac", "bc"), ifelse(act$kind == "alpha", "af", "bf"))
    act$idx <- NA_integer_
    for (pn in c("af", "bf")) { i <- act$par == pn; act$idx[i] <- seq_len(sum(i)) }
    for (pn in c("ac", "bc")) { i <- act$par == pn; act$idx[i] <- match(act$label[i], unique(act$label[i])) }
    act
  } else NULL
  S <- censor_build(S, sp$censor)                                       # censor = : limits and censored values
  S$par0 <- rtmb_start(S)
  S
}

rtmb_start <- function(S) {
  p <- list()
  pm <- pmin(pmax(colSums(S$Yb * S$Mb) / pmax(colSums(S$Mb), 1), 0.02), 0.98)
  p$d <- unname(stats::qlogis(pm))
  Yc <- S$Yc; Yc[S$Mc == 0] <- NA
  csd <- unname(apply(Yc, 2, stats::sd, na.rm = TRUE)); cm <- unname(colMeans(Yc, na.rm = TRUE))
  ld <- S$ld
  pos <- ld[ld$kind == "pos", ]; p$lpos <- if (nrow(pos)) {
    v <- tapply(ifelse(pos$type != "c", 0, log(0.5 * csd[ifelse(pos$type == "c", pos$j, 1)])), pos$idx, `[`, 1); unname(v) } else numeric(0)
  p$lreal <- rep(0, if (any(ld$kind == "real")) max(ld$idx[ld$kind == "real"]) else 0)
  p$xi <- cm; p$lsig <- log(0.8 * csd)
  p$eps <- rep(0, length(S$sh)); p$ldl <- rep(0, length(S$sh))
  p$kap <- rep(0, nrow(S$kap)); p$eta <- rep(0, nrow(S$eta))
  t1 <- numeric(); linc <- numeric(); pb <- numeric(); lr <- numeric()   # thresholds from the cumulative proportions
  for (j in seq_along(S$O)) {
    cp <- cumsum(colSums(S$ind[[j]]))[-S$K[j]] / sum(S$ind[[j]])
    bb <- stats::qlogis(pmin(pmax(cp, 0.01), 0.99)) / 1.2
    if (length(bb) > 1) for (k in 2:length(bb)) bb[k] <- max(bb[k], bb[k - 1] + 0.1)
    if (S$otype[j] == "grm") { t1 <- c(t1, bb[1]); if (length(bb) > 1) linc <- c(linc, log(diff(bb))) }
    else { pb <- c(pb, bb); if (S$otype[j] == "tppcm") lr <- c(lr, rep(0, S$K[j] - 2)) }
  }
  p$t1 <- t1; p$linc <- linc; p$pb <- pb; p$lr <- lr
  for (pn in c("af", "ac", "bf", "bc")) p[[pn]] <- rep(0, if (is.null(S$mod)) 0 else max(c(0, S$mod$idx[S$mod$par == pn])))
  for (i in seq_len(S$nF)) {
    fi <- S$f[[i]]; s <- paste0
    p[[s("beta", i)]] <- rep(0, length(fi$cols))
    if (fi$free_sd) {
      # factor with fixed loadings on continuous indicators: SD of the (scaled) person means
      r <- ld[ld$f == i, ]
      v <- if (all(r$kind == "fixed" & r$type == "c")) {
        Z <- sweep(sweep(Yc[, r$j, drop = FALSE], 2, cm[r$j]), 2, r$value, "/")
        max(stats::var(rowMeans(Z, na.rm = TRUE), na.rm = TRUE) - mean(csd[r$j]^2 / r$value^2) / nrow(r), 0.01)
      } else 0.25
      p[[s("lsd", i)]] <- 0.5 * log(v)
    }
    p[[s("alpha", i)]] <- rep(0, length(fi$scale))
    if (i == 2) {
      if (fi$path) p$gam <- 0
      if (fi$shash) { p$eps2 <- 0; p$ldl2 <- 0; p$eta2 <- rep(0, length(fi$skew)) }
      if (S$cor) p$atr <- 0
    }
  }
  p
}

# moderator sum: sum_k coef[k] * value of moderator k (a covariate column or a latent variable)
modsum <- function(coef, m, X, Fl, fs) {
  out <- 0
  for (k in seq_along(m)) out <- out + coef[k] * (if (m[k] %in% fs) Fl[[match(m[k], fs)]] else X[, m[k]])
  out
}

lin_x <- function(beta, cols, X) if (length(cols)) as.vector(X[, cols, drop = FALSE] %*% beta) else 0

# latent variables at standard-normal values z1, z2 for persons `rows`
lat_values <- function(p, S, z1, z2, rows) {
  X <- S$X[rows, , drop = FALSE]; f1i <- S$f[[1]]
  ls1 <- if (f1i$free_sd) p$lsd1 else log(f1i$sd)
  if (length(f1i$scale)) ls1 <- ls1 + modsum(p$alpha1, f1i$scale, X, NULL, S$fs)
  f1 <- lin_x(p$beta1, f1i$cols, X) + exp(ls1) * z1
  if (S$nF == 1) return(list(f1))
  f2i <- S$f[[2]]; Fl <- list(f1)
  mu2 <- lin_x(p$beta2, f2i$cols, X) + if (f2i$path) p$gam * f1 else 0
  ls2 <- if (f2i$free_sd) p$lsd2 else log(f2i$sd)
  if (length(f2i$scale)) ls2 <- ls2 + modsum(p$alpha2, f2i$scale, X, Fl, S$fs)
  u <- if (S$cor) { r <- tanh(p$atr); r * z1 + sqrt(1 - r * r) * z2 } else z2
  if (f2i$shash) {
    e2 <- p$eps2 + if (length(f2i$skew)) modsum(p$eta2, f2i$skew, X, Fl, S$fs) else 0
    u <- T_std(u, e2, exp(p$ldl2))
  }
  list(f1, mu2 + exp(ls2) * u)
}

loading <- function(p, r) switch(r$kind, fixed = r$value, pos = exp(p$lpos[r$idx]), real = p$lreal[r$idx])

# log-likelihood of the indicators, one value per row (person x node)
# moderation of indicator v: the loading multiplier exp(alpha'z) of its loading on f, and its shift delta'z
mod_eff <- function(p, S, X, v, kind, f = "") {
  if (is.null(S$mod)) return(NULL)
  m <- S$mod[S$mod$target == v & S$mod$kind == kind & (kind == "beta" | S$mod$factor == f), ]
  if (!nrow(m)) return(NULL)
  out <- 0
  for (i in seq_len(nrow(m))) out <- out + p[[m$par[i]]][m$idx[i]] * X[, m$mod[i]]
  out
}
# loading of row r at the rows of X (with its moderation)
loading_at <- function(p, S, r, X) {
  l <- loading(p, S$ld[r, ]); a <- mod_eff(p, S, X, S$ld$rhs[r], "alpha", S$ld$lhs[r])
  if (is.null(a)) l else l * exp(a)
}
# thresholds b_1 < b_2 < ... of ordinal indicator j
thresholds <- function(p, S, j) {
  ix <- S$tix[[j]]
  if (S$otype[j] != "grm") return(p$pb[ix$pb])
  b <- p$t1[ix$t1]
  if (S$K[j] > 2) for (k in 2:(S$K[j] - 1)) b <- c(b, b[k - 1] + exp(p$linc[ix$linc[k - 1]]))
  b
}
# step ratio r_k (k >= 2; r_1 = 1) of a tppcm item, 1 otherwise (AD-safe: no c() of numbers and advectors)
step_ratio <- function(p, S, j, k) if (S$otype[j] != "tppcm" || k == 1) 1 else exp(p$lr[S$tix[[j]]$lr[k - 1]])

ind_ll <- function(p, S, Fl, rows) {
  X <- S$X[rows, , drop = FALSE]; ld <- S$ld; ll <- 0
  for (j in seq_along(S$B)) {
    eta <- p$d[j]; L1 <- NULL
    for (r in S$ld_b[[j]]) { L <- loading_at(p, S, r, X); if (is.null(L1)) L1 <- L; eta <- eta + L * Fl[[ld$f[r]]] }
    dl <- mod_eff(p, S, X, S$B[j], "beta"); if (!is.null(dl)) eta <- eta - L1 * dl     # shift on the latent scale
    ll <- ll + S$Mb[rows, j] * RTMB::dbinom_robust(S$Yb[rows, j], 1, eta, log = TRUE)
  }
  for (j in seq_along(S$O)) {                                          # P(y > k) = logistic{a (f - delta - b_k)}
    r <- S$ld_o[[j]]; A <- loading_at(p, S, r, X); th <- Fl[[ld$f[r]]]
    dl <- mod_eff(p, S, X, S$O[j], "beta"); if (!is.null(dl)) th <- th - dl
    b <- thresholds(p, S, j); I <- S$ind[[j]][rows, , drop = FALSE]
    if (S$otype[j] != "grm") {                                         # adjacent-category logits: log weight lp_k, normalizer acc
      lp <- 0 * th; num <- 0 * th; acc <- 0 * th
      for (k in 2:S$K[j]) {
        lp <- lp + A * step_ratio(p, S, j, k - 1) * (th - b[k - 1]); num <- num + I[, k] * lp
        acc <- acc + log1p(exp(lp - acc))
      }
      ll <- ll + S$Mo[rows, j] * (num - acc)
      next
    }
    up <- 1; pr <- 0
    for (k in seq_len(S$K[j])) {
      lo <- if (k < S$K[j]) 1 / (1 + exp(-A * (th - b[k]))) else 0
      pr <- pr + I[, k] * (up - lo); up <- lo
    }
    ll <- ll + S$Mo[rows, j] * log(pr + (1 - S$Mo[rows, j]))
  }
  for (j in seq_along(S$C)) {
    mu <- p$xi[j]
    for (r in S$ld_c[[j]]) mu <- mu + loading_at(p, S, r, X) * Fl[[ld$f[r]]]
    dl <- mod_eff(p, S, X, S$C[j], "beta"); if (!is.null(dl)) mu <- mu + dl
    lsc <- p$lsig[j]; kr <- which(S$kap$j == j)
    if (length(kr)) lsc <- lsc + modsum(p$kap[kr], S$kap$m[kr], X, Fl, S$fs)
    x <- (S$Yc[rows, j] - mu) / exp(lsc)
    s <- match(j, S$sh)
    e <- NULL
    ld_x <- if (is.na(s)) -0.5 * x * x - 0.9189385332046727 else {
      er <- which(S$eta$j == j)
      e <- p$eps[s] + if (length(er)) modsum(p$eta[er], S$eta$m[er], X, Fl, S$fs) else 0
      ld_std(x, e, exp(p$ldl[s]))
    }
    cr <- if (is.null(S$cens)) integer() else which(S$cens$j == j)
    if (!length(cr)) { ll <- ll + S$Mc[rows, j] * (ld_x - lsc); next }
    # censored values: log F(lower) / log S(upper) instead of the log density (R/censor.R)
    L <- S$Lc[rows, j]; R <- S$Rc[rows, j]; dl <- if (is.na(s)) 1 else exp(p$ldl[s])
    lj <- (1 - L - R) * (ld_x - lsc)
    if (is.finite(S$cens$lower[cr])) lj <- lj + L * cens_logp(S$cens$lower[cr], mu, lsc, e, dl, upper = FALSE)
    if (is.finite(S$cens$upper[cr])) lj <- lj + R * cens_logp(S$cens$upper[cr], mu, lsc, e, dl, upper = TRUE)
    ll <- ll + S$Mc[rows, j] * lj
  }
  ll
}

# Laplace: joint negative log density with z1, z2 as random effects
make_joint <- function(S) {
  rows <- seq_len(S$N)
  function(p) {
    z2 <- if (S$nF == 2) p$z2 else 0
    llp <- ind_ll(p, S, lat_values(p, S, p$z1, z2, rows), rows) - 0.5 * p$z1^2 - 0.9189385332046727
    if (S$nF == 2) llp <- llp - 0.5 * p$z2^2 - 0.9189385332046727
    RTMB::REPORT(llp)
    -sum(llp) - if (is.null(S$prior)) 0 else S$prior$lp(p)
  }
}

# AGHQ: marginal negative log-likelihood with the nodes ND (data)
make_aghq <- function(S, ND) {
  N <- S$N; G <- ND$G
  function(p) {
    Fl <- lat_values(p, S, ND$z1, ND$z2, ND$rows)
    ll <- ind_ll(p, S, Fl, ND$rows) + ND$lw
    lacc <- ll[seq_len(N)]                           # per-person log-sum-exp over the nodes
    for (h in seq_len(G)[-1]) lacc <- RTMB::logspace_add(lacc, ll[(h - 1) * N + seq_len(N)])
    F1 <- Fl[[1]]; RTMB::REPORT(ll); RTMB::REPORT(F1)
    if (S$nF == 2) { F2 <- Fl[[2]]; RTMB::REPORT(F2) }
    -sum(lacc + ND$cst) - if (is.null(S$prior)) 0 else S$prior$lp(p)
  }
}

re_names <- function(S) if (S$nF == 2) c("z1", "z2") else "z1"

# conditional modes and Hessians of (z1, z2) at the fixed parameters pf
centers <- function(S, pf, z0) {
  obj <- RTMB::MakeADFun(make_joint(S), c(pf, z0), random = re_names(S), silent = TRUE)
  obj$fn(obj$par)
  centers_at(S, obj)
}
# ... from a Laplace object just evaluated at the fixed parameters of interest
centers_at <- function(S, obj) {
  lp <- obj$env$last.par; zr <- lp[obj$env$random]; N <- S$N
  H <- obj$env$spHess(lp, random = TRUE)
  cst <- obj$report(lp)$llp
  i <- seq_len(N)
  # guard against persons whose inner problem is not (yet) well behaved
  if (S$nF == 1) return(list(m1 = zr, l11 = 1 / sqrt(pmax(Matrix::diag(H), 0.01)), cst = cst, z = list(z1 = zr)))
  a <- pmax(H[cbind(i, i)], 0.01); b <- H[cbind(i, N + i)]; c <- pmax(H[cbind(N + i, N + i)], 0.01)
  b <- ifelse(a * c - b * b > 1e-4 * a * c, b, 0)
  det <- a * c - b * b
  s11 <- c / det; s21 <- -b / det; s22 <- a / det                  # inverse Hessian
  l11 <- sqrt(s11); l21 <- s21 / l11; l22 <- sqrt(pmax(s22 - l21^2, 1e-12))
  list(m1 = zr[i], m2 = zr[N + i], l11 = l11, l21 = l21, l22 = l22, cst = cst, z = list(z1 = zr[i], z2 = zr[N + i]))
}

node_data <- function(S, cn, k) {
  gq <- statmod::gauss.quad(k, "hermite"); N <- S$N; s2 <- sqrt(2)
  if (S$nF == 1) {
    G <- k; g1 <- rep(gq$nodes, each = N); lwg <- rep(log(gq$weights) + gq$nodes^2, each = N)
    z1 <- rep(cn$m1, G) + s2 * rep(cn$l11, G) * g1
    lw <- lwg + rep(log(cn$l11) + 0.5 * log(2), G) - 0.5 * z1^2 - 0.9189385332046727 - rep(cn$cst, G)
    return(list(G = G, rows = rep(seq_len(N), G), z1 = z1, z2 = 0, lw = lw, cst = cn$cst))
  }
  gg <- as.matrix(expand.grid(seq_len(k), seq_len(k))); G <- nrow(gg)
  g1 <- rep(gq$nodes[gg[, 1]], each = N); g2 <- rep(gq$nodes[gg[, 2]], each = N)
  lwg <- rep(log(gq$weights[gg[, 1]]) + log(gq$weights[gg[, 2]]) + gq$nodes[gg[, 1]]^2 + gq$nodes[gg[, 2]]^2, each = N)
  z1 <- rep(cn$m1, G) + s2 * rep(cn$l11, G) * g1
  z2 <- rep(cn$m2, G) + s2 * (rep(cn$l21, G) * g1 + rep(cn$l22, G) * g2)
  lw <- lwg + rep(log(cn$l11) + log(cn$l22) + log(2), G) - 0.5 * (z1^2 + z2^2) - 2 * 0.9189385332046727 - rep(cn$cst, G)
  list(G = G, rows = rep(seq_len(N), G), z1 = z1, z2 = z2, lw = lw, cst = cn$cst)
}

# SHASH shapes bounded: tail weight delta in [1/10, 10], |skewness eps| <= 5, |moderation| <= 5
rtmb_bounds <- function(par, side) {
  n <- names(par); b <- rep(if (side == "lo") -Inf else Inf, length(par))
  lim <- c(ldl = log(10), ldl2 = log(10), eps = 5, eps2 = 5, eta = 5, eta2 = 5)
  i <- n %in% names(lim); b[i] <- (if (side == "lo") -1 else 1) * lim[n[i]]
  b
}

fam_pars <- c("eps", "ldl", "kap", "eta", "alpha1", "alpha2", "eps2", "ldl2", "eta2")

fit_rtmb <- function(S, ctrl, progress = TRUE) {
  t0 <- Sys.time()
  say <- function(...) if (progress) message(...)
  meth <- if (ctrl$method == "elgm") ctrl$inner else ctrl$method      # elgm: the posterior mode first
  bnd <- rtmb_bounds
  nlm <- function(obj) withCallingHandlers(   # NaN trial points in the line search are rejected by nlminb
    stats::nlminb(obj$par, obj$fn, obj$gr, lower = bnd(obj$par, "lo"), upper = bnd(obj$par, "up"),
                  control = list(eval.max = 3000, iter.max = 2000, trace = as.integer(ctrl$trace))),
    warning = function(w) if (grepl("NA/NaN function evaluation", conditionMessage(w))) invokeRestart("muffleWarning"))
  N <- S$N; pf <- S$par0
  z0 <- list(z1 = rep(0, N)); if (S$nF == 2) z0$z2 <- rep(0, N)
  # stage 1: normal, homoscedastic; stage 2: the full model (Laplace); then AGHQ from there
  fp <- intersect(fam_pars, names(pf)); fp <- fp[vapply(pf[fp], length, 1L) > 0]
  stages <- if (length(fp)) list(lapply(pf[fp], function(v) factor(rep(NA, length(v)))), list()) else list(list())
  for (st in seq_along(stages)) {
    say(sprintf("Laplace%s ...", if (length(stages) > 1) c(" (normal start)", "")[st] else ""))
    objL1 <- RTMB::MakeADFun(make_joint(S), c(pf, z0), random = re_names(S), map = stages[[st]], silent = TRUE)
    opt1 <- tryCatch(nlm(objL1), error = function(e) NULL)
    if (is.null(opt1) || !is.finite(opt1$objective)) {
      if (st == 1) stop("the Laplace start failed; check the model and the data", call. = FALSE)
      if (meth == "laplace") stop("the Laplace fit of the full model failed; try method = \"aghq\"", call. = FALSE)
      say("Laplace fit of the full model failed; AGHQ starts from the normal model")
      break
    }
    objL <- objL1; opt <- opt1
    pl <- objL$env$parList(opt$par); pf <- pl[names(S$par0)]; z0 <- pl[re_names(S)]
  }
  res <- list(method = meth, laplace = list(objective = opt$objective, convergence = opt$convergence))
  k <- ctrl$k %||% if (S$nF == 1) 21L else 9L
  if (meth == "laplace") {
    sdr <- RTMB::sdreport(objL)
    V <- sdr$cov.fixed; obj_value <- opt$objective; conv <- opt$convergence
    cn <- centers(S, pf, z0); ND <- node_data(S, cn, k)
    objA <- RTMB::MakeADFun(make_aghq(S, ND), pf, silent = TRUE)
    rounds <- opt$objective; xA <- objA$par
    gr <- objL$gr(opt$par); ab <- rep(FALSE, length(xA))
  } else {
    # AGHQ with the nodes at the conditional modes; after each optimization the nodes are
    # recentred at the new estimates and the objective re-evaluated there (the fixed-node
    # objective can be exploited by large steps); a worse value halves the step.
    at <- function(pf, z0) {
      cn <- centers(S, pf, z0); ND <- node_data(S, cn, k)
      obj <- RTMB::MakeADFun(make_aghq(S, ND), pf, silent = TRUE)
      list(pf = pf, cn = cn, ND = ND, obj = obj, value = obj$fn(obj$par))
    }
    mix <- function(a, b, t) Map(function(u, v) u + t * (v - u), a, b)
    cur <- at(pf, z0); rounds <- cur$value; conv <- 0
    say(sprintf("AGHQ (k = %d) at the Laplace estimates: -2logL %.3f", k, 2 * cur$value))
    for (r in seq_len(ctrl$max_rounds)) {
      optA <- nlm(cur$obj); conv <- optA$convergence
      pnew <- cur$obj$env$parList(optA$par)
      new <- tryCatch(at(pnew, cur$cn$z), error = function(e) NULL)
      if (is.null(new) || !is.finite(new$value) || new$value > cur$value + ctrl$tol / 2) {
        new <- NULL
        for (t in c(0.5, 0.25, 0.125, 0.0625)) {
          tr <- tryCatch(at(mix(cur$pf, pnew, t), cur$cn$z), error = function(e) NULL)
          if (!is.null(tr) && is.finite(tr$value) && tr$value < cur$value) { new <- tr; break }
        }
        if (is.null(new)) {
          warning("AGHQ: no improvement over the previous round; the model may be poorly identified for these data", call. = FALSE)
          break
        }
      }
      done <- abs(cur$value - new$value) < ctrl$tol / 2
      cur <- new; rounds <- c(rounds, cur$value)
      say(sprintf("AGHQ (k = %d) round %d: -2logL %.3f", k, r, 2 * cur$value))
      if (done) break
    }
    pf <- cur$pf; ND <- cur$ND; objA <- cur$obj
    xA <- objA$par; obj_value <- cur$value
    # polish at the final nodes (the recentred objective has a small gradient at the old optimum)
    optP <- tryCatch(nlm(objA), error = function(e) NULL)
    if (!is.null(optP) && is.finite(optP$objective) && optP$objective <= obj_value + ctrl$tol) {
      xA <- optP$par; pf <- objA$env$parList(xA)
    }
    gr_conv <- objA$gr(xA)                       # convergence: gradient of the objective that was optimized
    # final value with the nodes recentred at the estimates, and a check with k + 6 nodes
    cnf <- centers(S, pf, cur$cn$z); ND <- node_data(S, cnf, k)
    objA <- RTMB::MakeADFun(make_aghq(S, ND), pf, silent = TRUE); xA <- objA$par; obj_value <- objA$fn(xA)
    chk <- RTMB::MakeADFun(make_aghq(S, node_data(S, cnf, k + 6)), pf, silent = TRUE)
    quad <- c(obj_value, chk$fn(chk$par)); names(quad) <- c(sprintf("k=%d", k), sprintf("k=%d", k + 6))
    if (abs(2 * diff(quad)) > 0.5)
      warning(sprintf("AGHQ: -2logL changes by %.2f with %d instead of %d nodes per dimension; increase k (rtmb_control(k = %d))",
                      2 * diff(quad), k + 6, k, k + 6), call. = FALSE)
    H <- tryCatch(objA$he(xA), error = function(e) NULL)
    if (is.null(H) || any(!is.finite(H))) H <- stats::optimHess(xA, objA$fn, objA$gr)
    H <- (H + t(H)) / 2
    gr <- gr_conv
    # parameters at a bound are treated as fixed (SEs conditional on them)
    ab <- abs(xA - bnd(xA, "lo")) < 1e-6 | abs(xA - bnd(xA, "up")) < 1e-6
    if (any(ab)) {
      warning("SHASH shape parameters at their bounds (", paste(unique(names(xA)[ab]), collapse = ", "),
              "): this residual is close to normal or the tail weight is not identified; the standard errors condition on them", call. = FALSE)
      gr[ab] <- 0
    }
    V <- matrix(0, length(xA), length(xA))
    Vf <- tryCatch(solve(H[!ab, !ab, drop = FALSE]), error = function(e) NULL)
    if (is.null(Vf) || any(!is.finite(Vf)) || any(diag(Vf) < 0)) {
      warning("the Hessian is not positive definite: standard errors are not available for all parameters", call. = FALSE)
      ev <- eigen(H[!ab, !ab, drop = FALSE], symmetric = TRUE)
      Vf <- ev$vectors %*% diag(ifelse(ev$values > 1e-8, 1 / ev$values, NA), length(ev$values)) %*% t(ev$vectors)
    }
    V[!ab, !ab] <- Vf
    big <- which(!ab & is.finite(diag(V)) & sqrt(pmax(diag(V), 0)) > 50)
    if (length(big))
      warning("very large standard errors (internal scale) for ", paste(unique(names(xA)[big]), collapse = ", "),
              ": the model is close to non-identified for these data, or a parameter diverges (e.g. a discrimination or an intercept running off)", call. = FALSE)
  }
  rep <- objA$report(xA)
  W <- matrix(rep$ll, N); W <- exp(W - apply(W, 1, max)); W <- W / rowSums(W)
  sc <- lapply(seq_len(S$nF), function(i) {
    Fm <- matrix(rep[[paste0("F", i)]], N, ND$G); m <- rowSums(W * Fm)
    list(mean = m, sd = sqrt(pmax(rowSums(W * Fm^2) - m^2, 0)))
  })
  c(res, list(par = pf, x = xA, vcov = V, objective = obj_value, convergence = conv, rounds = rounds, k = k, at_bound = names(xA)[ab],
              centers = if (meth == "aghq") cnf else cn,                    # node placement of the final objective (scores)
              quad_check = if (meth == "aghq") 2 * quad, ab = ab,
              maxgrad = max(abs(gr)), npar = length(xA), scores = sc,
              secs = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

# ---- reported quantities (natural scale) ---------------------------------------------------
# Returns a function of the parameter list giving a named vector; with numeric parameters it
# gives the estimates and their names, under RTMB the same function is differentiated for
# delta-method standard errors.
natural_fun <- function(S) {
  function(p) {
    out <- list(); add <- function(n, v) out[[n]] <<- v
    ld <- S$ld
    for (r in seq_len(nrow(ld))) add(sprintf("%s|=~|%s", ld$lhs[r], ld$rhs[r]), loading(p, ld[r, ]))
    for (i in seq_len(S$nF)) {
      fi <- S$f[[i]]; f <- fi$name
      if (i == 2 && fi$path) add(sprintf("%s|~|%s", f, S$fs[1]), p$gam)
      for (k in seq_along(fi$covs)) add(sprintf("%s|~|%s", f, fi$covs[k]), p[[paste0("beta", i)]][k])
    }
    if (S$cor) add(sprintf("%s|~~|%s", S$fs[1], S$fs[2]), tanh(p$atr))
    for (i in seq_len(S$nF)) {
      fi <- S$f[[i]]; f <- fi$name
      ls <- if (fi$free_sd) p[[paste0("lsd", i)]] else log(fi$sd)
      if (!fi$shash && !length(fi$scale)) add(sprintf("%s|~~|%s", f, f), exp(2 * ls))
      else {
        add(sprintf("%s|sd|", f), exp(ls))
        for (k in seq_along(fi$scale)) add(sprintf("V(%s)|~|%s", f, fi$scale[k]), p[[paste0("alpha", i)]][k])
        if (fi$shash) {
          add(sprintf("%s|skew|", f), p$eps2)
          for (k in seq_along(fi$skew)) add(sprintf("%s|skew~|%s", f, fi$skew[k]), p$eta2[k])
          add(sprintf("%s|tail|", f), exp(p$ldl2))
        }
      }
    }
    for (j in seq_along(S$B)) add(sprintf("%s|~1|", S$B[j]), p$d[j])
    for (j in seq_along(S$O)) {
      b <- thresholds(p, S, j); for (k in seq_along(b)) add(sprintf("%s|thr|t%d", S$O[j], k), b[k])
      if (S$otype[j] == "tppcm" && S$K[j] > 2) { a <- loading(p, ld[S$ld_o[[j]], ])
        for (k in 2:(S$K[j] - 1)) add(sprintf("%s|stepa|a%d", S$O[j], k), a * step_ratio(p, S, j, k)) }
    }
    if (!is.null(S$mod)) for (i in seq_len(nrow(S$mod))) { m <- S$mod[i, ]
      add(sprintf("E(%s)|~|%s", m$target, if (m$kind == "alpha") paste0(m$mod, ":", m$factor) else m$mod), p[[m$par]][m$idx]) }
    for (j in seq_along(S$C)) add(sprintf("%s|~1|", S$C[j]), p$xi[j])
    for (j in seq_along(S$C)) {
      v <- S$C[j]; kr <- which(S$kap$j == j); s <- match(j, S$sh)
      if (is.na(s) && !length(kr)) { add(sprintf("%s|~~|%s", v, v), exp(2 * p$lsig[j])); next }
      add(sprintf("%s|sd|", v), exp(p$lsig[j]))
      for (r in kr) add(sprintf("V(%s)|~|%s", v, S$kap$m[r]), p$kap[r])
      if (!is.na(s)) {
        add(sprintf("%s|skew|", v), p$eps[s])
        for (r in which(S$eta$j == j)) add(sprintf("%s|skew~|%s", v, S$eta$m[r]), p$eta[r])
        add(sprintf("%s|tail|", v), exp(p$ldl[s]))
      }
    }
    # IRT parameterization of binary items with one loading: a = lambda, b = -d / lambda
    for (j in seq_along(S$B)) if (length(S$ld_b[[j]]) == 1) {
      a <- loading(p, ld[S$ld_b[[j]], ]); add(sprintf("%s|irt_a|", S$B[j]), a); add(sprintf("%s|irt_b|", S$B[j]), -p$d[j] / a)
    }
    out
  }
}

# Estimates and SEs of fun(par) (a named list of scalars) at the fit: the delta method (maximum
# likelihood), or the posterior mean, SD and 95% interval under the outer quadrature (ELGM).
# tr: scale of the ELGM interval per output ("log", "atanh" or "id").
# delta method (ML) or ELGM posterior summaries of functions of the parameters; cov = TRUE adds
# their joint covariance: J V J' (ML), or the covariance over the ELGM nodes plus J Sc J' (the
# same decomposition as the posterior SDs of elgm_summary())
delta <- function(fit, fun, tr = NULL, cov = FALSE) {
  pl <- fit$rtmb$par
  lst <- fun(pl); est <- vapply(lst, function(e) as.numeric(e)[1], 0)
  obj <- RTMB::MakeADFun(function(p) { v <- do.call(c, lapply(fun(p), RTMB::AD)); RTMB::ADREPORT(v); 0 },
                         pl, ADreport = TRUE, silent = TRUE)
  J <- obj$gr(fit$rtmb$x); if (is.null(dim(J))) J <- matrix(J, nrow = length(est))
  fixed <- rowSums(abs(J)) == 0
  E <- fit$rtmb$elgm
  if (is.null(E)) {
    se <- sqrt(pmax(rowSums((J %*% fit$rtmb$vcov) * J), 0)); se[fixed] <- NA
    out <- list(est = unname(est), se = se, lo = est - 1.959964 * se, hi = est + 1.959964 * se, names = names(lst))
    if (cov) { out$V <- J %*% fit$rtmb$vcov %*% t(J); out$J <- J }
    return(out)
  }
  G <- vapply(seq_len(nrow(E$nodes)), function(g) {
    if (E$prob[g] == 0) return(est)
    vapply(fun(vec2list(E$nodes[g, ], pl)), function(e) as.numeric(e)[1], 0)
  }, est)
  if (is.null(dim(G))) G <- matrix(G, nrow = length(est))
  s <- elgm_summary(E, G, J, tr %||% rep("id", length(est)))
  s$se[fixed] <- NA; s$lo[fixed] <- s$hi[fixed] <- NA
  if (cov) { Gc <- G - as.vector(G %*% E$prob); s$V <- Gc %*% (E$prob * t(Gc)) + J %*% E$Sc %*% t(J) }
  c(s, list(names = names(lst)))
}

rtmb_estimates <- function(fit) {
  S <- fit$S; nf <- natural_fun(S)
  nm <- names(nf(fit$rtmb$par))
  parts <- do.call(rbind, strsplit(sub("\\|$", "| ", nm), "|", fixed = TRUE))
  op <- parts[, 2]; lhs <- parts[, 1]; rhs <- trimws(parts[, 3]); op[op == "thr"] <- "|"; op[op == "stepa"] <- "|a"
  pos_ld <- paste(S$ld$lhs, S$ld$rhs)[S$ld$kind == "pos"]
  tr <- ifelse(op %in% c("sd", "tail") | (op == "~~" & lhs == rhs) | (op == "=~" & paste(lhs, rhs) %in% pos_ld), "log",
               ifelse(op == "~~", "atanh", "id"))
  d <- delta(fit, nf, tr)
  e <- data.frame(lhs = lhs, op = op, rhs = rhs, est = d$est, se = d$se, stringsAsFactors = FALSE)
  bayes <- !is.null(fit$rtmb$elgm)
  e$z <- if (bayes) NA_real_ else e$est / e$se
  e$pvalue <- if (bayes) NA_real_ else 2 * stats::pnorm(-abs(e$z))
  e$ci.lower <- d$lo; e$ci.upper <- d$hi
  e$fixed <- is.na(e$se)
  if (bayes) { e$prior <- unname(S$prior$txt[nm]); e$prior[is.na(e$prior)] <- "" }
  e
}

# ---- quantile effects ----------------------------------------------------------------------
# Population rows for the expectation over the moderators: persons (when covariates matter)
# x Gauss-Hermite nodes (probabilists') of the needed standard normals.
pop_rows <- function(S, need_x, need_z1, need_z2, g1 = 15, g2 = 9, only = NULL) {
  rows <- if (need_x) (only %||% seq_len(S$N)) else 1L
  gq1 <- if (need_z2) statmod::gauss.quad.prob(g2, "normal") else statmod::gauss.quad.prob(g1, "normal")
  z1 <- if (need_z1) gq1$nodes else 0; w1 <- if (need_z1) gq1$weights else 1
  z2 <- if (need_z2) gq1$nodes else 0; w2 <- if (need_z2) gq1$weights else 1
  g <- expand.grid(r = seq_along(rows), a = seq_along(z1), b = seq_along(z2))
  list(rows = rows[g$r], z1 = z1[g$a], z2 = z2[g$b], w = w1[g$a] * w2[g$b] / length(rows))
}

# standardized SHASH transform T_std(u; e, dl): value, d/du, d/de
shash_T <- function(u, e, dl) {
  a <- eff_a(e, dl); da <- 1 / (dl * cosh(e / dl / 8)^2)   # d a / d e
  w <- asinh(u) / dl + a; k <- cosh(a) / dl
  Tv <- (sinh(w) - sinh(a)) / k
  list(T = Tv, Tu = cosh(w) / (dl * sqrt(1 + u * u)) / k,
       Te = da * ((cosh(w) - cosh(a)) / k - Tv * tanh(a)))
}

quantile_targets <- function(S) {
  out <- list()
  if (S$nF == 2) {
    fi <- S$f[[2]]; f1 <- S$fs[1]
    pr <- unique(c(if (fi$path || S$cor || f1 %in% c(fi$scale, fi$skew)) f1, fi$covs, setdiff(c(fi$scale, fi$skew), f1)))
    if (length(pr)) out[[fi$name]] <- list(type = "latent", pred = pr)
  }
  for (j in seq_along(S$C)) {
    m <- unique(c(S$kap$m[S$kap$j == j], S$eta$m[S$eta$j == j]))
    em <- if (is.null(S$mod)) character() else unique(S$mod$mod[S$mod$target == S$C[j]])      # E(y) ~ z, E(y) ~ z:f
    af <- if (is.null(S$mod)) character() else unique(S$mod$factor[S$mod$target == S$C[j] & S$mod$kind == "alpha"])
    pr <- unique(c(S$ld$lhs[S$ld_c[[j]]], m, em))
    out[[S$C[j]]] <- list(type = "indicator", j = j, pred = pr, mods = unique(c(m, em, af)))
  }
  out
}

qeffect_fun <- function(S, P, only = NULL) {           # only: average over these persons (default: all)
  zp <- stats::qnorm(P); fs <- S$fs; tg <- quantile_targets(S)
  # population for each target
  pops <- lapply(tg, function(t) {
    if (t$type == "latent") {
      fi <- S$f[[2]]
      nonlin <- fi$shash || length(fi$scale) > 0 || S$cor
      pop_rows(S, need_x = nonlin && length(S$Xn) > 0, need_z1 = nonlin && (S$cor || fs[1] %in% c(fi$scale, fi$skew)), need_z2 = FALSE, only = only)
    } else {
      fm <- intersect(t$mods, fs)
      pop_rows(S, need_x = length(t$mods) > 0 && length(S$Xn) > 0, need_z1 = length(fm) > 0,
               need_z2 = length(fs) == 2 && fs[2] %in% fm, only = only)
    }
  })
  f <- function(p) {
    out <- list()
    for (n in names(tg)) {
      t <- tg[[n]]; po <- pops[[n]]; X <- S$X[po$rows, , drop = FALSE]
      Fl <- lat_values(p, S, po$z1, po$z2, po$rows)
      res <- vector("list", length(P))
      for (q in seq_along(P)) {
        if (t$type == "latent") {
          fi <- S$f[[2]]; f1i <- S$f[[1]]
          ls1 <- if (f1i$free_sd) p$lsd1 else log(f1i$sd)
          if (length(f1i$scale)) ls1 <- ls1 + modsum(p$alpha1, f1i$scale, X, NULL, fs)
          z1 <- po$z1; s1 <- exp(ls1)
          ls2 <- if (fi$free_sd) p$lsd2 else log(fi$sd)
          if (length(fi$scale)) ls2 <- ls2 + modsum(p$alpha2, fi$scale, X, Fl, fs)
          r <- if (S$cor) tanh(p$atr) else 0
          u <- r * z1 + sqrt(1 - r * r) * zp[q]
          if (fi$shash) {
            e2 <- p$eps2 + if (length(fi$skew)) modsum(p$eta2, fi$skew, X, Fl, fs) else 0
            Tt <- shash_T(u, e2, exp(p$ldl2))
          } else Tt <- list(T = u, Tu = 1, Te = 0)
          s2 <- exp(ls2); dl <- list()
          for (k in t$pred) {
            g <- 0
            if (k == fs[1]) {
              if (fi$path) g <- g + p$gam
              g <- g + s2 * Tt$Tu * r / s1                                  # through z1 = (f1 - mu1) / s1
            } else {
              b <- match(k, fi$covs); if (!is.na(b)) g <- g + p$beta2[b]
              b1 <- match(k, S$f[[1]]$covs); a1 <- match(k, f1i$scale)       # z1 moves with x for fixed f1
              dz1 <- (if (!is.na(b1)) -p$beta1[b1] / s1 else 0) + (if (!is.na(a1)) -z1 * p$alpha1[a1] else 0)
              g <- g + s2 * Tt$Tu * r * dz1
            }
            a <- match(k, fi$scale); if (!is.na(a)) g <- g + s2 * p$alpha2[a] * Tt$T
            e <- match(k, fi$skew); if (!is.na(e)) g <- g + s2 * p$eta2[e] * Tt$Te
            dl[[k]] <- sum(po$w * g)
          }
          res[[q]] <- dl
        } else {
          j <- t$j; s <- match(j, S$sh); kr <- which(S$kap$j == j); er <- which(S$eta$j == j)
          lsc <- p$lsig[j] + if (length(kr)) modsum(p$kap[kr], S$kap$m[kr], X, Fl, fs) else 0
          if (is.na(s)) Tt <- list(T = zp[q], Te = 0) else {
            e <- p$eps[s] + if (length(er)) modsum(p$eta[er], S$eta$m[er], X, Fl, fs) else 0
            Tt <- shash_T(zp[q], e, exp(p$ldl[s]))
          }
          sc <- exp(lsc); dl <- list()
          for (k in t$pred) {
            g <- 0
            for (r in S$ld_c[[j]]) if (S$ld$lhs[r] == k) g <- g + sum(po$w * loading_at(p, S, r, X))   # average loading (moderated)
            if (!is.null(S$mod)) {                                                                    # shifts and loading moderation by k
              mb <- S$mod[S$mod$target == S$C[j] & S$mod$kind == "beta" & S$mod$mod == k, ]
              for (i in seq_len(nrow(mb))) g <- g + p[[mb$par[i]]][mb$idx[i]]
              ma <- S$mod[S$mod$target == S$C[j] & S$mod$kind == "alpha" & S$mod$mod == k, ]
              for (i in seq_len(nrow(ma))) { r <- S$ld_c[[j]][S$ld$lhs[S$ld_c[[j]]] == ma$factor[i]]
                g <- g + sum(po$w * p[[ma$par[i]]][ma$idx[i]] * loading_at(p, S, r, X) * Fl[[match(ma$factor[i], fs)]]) }
            }
            a <- match(k, S$kap$m[kr]); if (!is.na(a)) g <- g + sum(po$w * sc * p$kap[kr[a]] * Tt$T)
            b <- match(k, S$eta$m[er]); if (!is.na(b)) g <- g + sum(po$w * sc * p$eta[er[b]] * Tt$Te)
            dl[[k]] <- g
          }
          res[[q]] <- dl
        }
      }
      for (q in seq_along(P)) for (k in t$pred) out[[sprintf("%s|%s|%s", n, k, format(P[q]))]] <- res[[q]][[k]]
    }
    out
  }
  list(f = f, n = sum(vapply(tg, function(t) length(t$pred), 1L)) * length(P))
}

#' Quantile effects of an RTMB fit
#'
#' For a model fitted with `engine = "rtmb"`: the effect of each predictor on the conditional
#' p-quantile of the second latent variable (e.g. `speed ~ theta + x`) and of each continuous
#' indicator (e.g. the cross-relation of `theta` with a log response time), averaged over the
#' model-implied distribution of the moderators (average marginal effects), with
#' delta-method standard errors (with `method = "elgm"`: posterior means, SDs and 95\% credible
#' intervals). With normal residuals and no moderators the effects are the
#' same at every level (the mean-model coefficients); `shash()` and `scale` / `skew`
#' moderators let them change with p. These are the counterparts of the coefficients of
#' `ald(p)` fits (engine `"jags"`), from one fit for all p.
#'
#' `band.lower` / `band.upper` form a simultaneous 95\% band over the levels `p` of each
#' target-predictor curve (sup-t band; Montiel Olea & Plagborg-Moller, 2019): every effect of the
#' curve lies in its band with probability 0.95, so the levels where the band excludes 0 can be
#' read off together. The critical value is the 0.95 quantile of max |Z| for Z normal with the
#' correlation of the curve's estimates (the joint delta-method covariance; with ELGM, the
#' posterior covariance). [quantile_test()] tests whether a curve changes with p.
#' @param object A `birt_rtmb` fit.
#' @param p Quantile levels.
#' @param ... Unused.
#' @return A data frame with `target`, `predictor`, `p`, `est`, `se`, `ci.lower`, `ci.upper`
#'   (pointwise), `band.lower`, `band.upper` (simultaneous over `p`); attribute `vcov`: the joint
#'   covariance of `est`.
#' @references Montiel Olea, J. L., & Plagborg-Moller, M. (2019). Simultaneous confidence bands:
#'   Theory, implementation, and an application to SVARs. *Journal of Applied Econometrics, 34*(1),
#'   1-17. \doi{10.1002/jae.2656}
#' @export
quantiles <- function(object, ...) UseMethod("quantiles")

#' @rdname quantiles
#' @export
quantiles.birt_rtmb <- function(object, p = object$control$p, ...) {
  qf <- qeffect_fun(object$S, p)
  if (!qf$n) return(data.frame(target = character(), predictor = character(), p = numeric(), est = numeric(), se = numeric()))
  d <- delta(object, qf$f, cov = TRUE)
  parts <- do.call(rbind, strsplit(d$names, "|", fixed = TRUE))
  out <- data.frame(target = parts[, 1], predictor = parts[, 2], p = as.numeric(parts[, 3]), est = d$est, se = d$se,
                    ci.lower = d$lo, ci.upper = d$hi, band.lower = NA_real_, band.upper = NA_real_)
  for (g in split(seq_len(nrow(out)), paste(out$target, out$predictor))) {
    ok <- g[!is.na(out$se[g]) & out$se[g] > 0]
    if (!length(ok)) next
    cv <- supt_crit(d$V[ok, ok, drop = FALSE])
    out$band.lower[ok] <- out$est[ok] - cv * out$se[ok]; out$band.upper[ok] <- out$est[ok] + cv * out$se[ok]
  }
  attr(out, "vcov") <- structure(d$V, dimnames = list(d$names, d$names))
  out
}

# 0.95 quantile of max_k |Z_k|, Z ~ N(0, R) with R the correlation of V (sup-t critical value;
# Montiel Olea & Plagborg-Moller, 2019); simulated with a fixed seed, 1.96 for one level
supt_crit <- function(V, level = 0.95, n = 1e5) {
  if (nrow(V) == 1) return(stats::qnorm(1 - (1 - level) / 2))
  s <- sqrt(diag(V)); R <- V / outer(s, s)
  e <- eigen(R, symmetric = TRUE); L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow(R))
  local_seed(20261005)
  Z <- matrix(stats::rnorm(n * nrow(R)), n) %*% t(L)
  unname(stats::quantile(apply(abs(Z), 1, max), level))
}

#' Does a quantile effect change with the quantile level?
#'
#' For a `birt_rtmb` fit: for each target and predictor of [quantiles()] (e.g. the cross-relation
#' `t ~ speed + ability` of a log response time, or `speed ~ ability + x`), a test of whether the
#' effect changes with the quantile level p, and the difference of each level from the median.
#'
#' In the model the location of the target is linear in the predictor and its scale and skewness
#' may depend on it, so the p-quantile effect changes with p exactly when the predictor moderates
#' the scale (`V(t) ~ k`) or the skewness of the target; equal effects at every level is the null
#' hypothesis of no heteroscedasticity, as in the test of equal slopes across regression quantiles
#' (Koenker & Bassett, 1982). The test is therefore the Wald test of these moderation parameters
#' (df = their number). A test of the differences of the effects across p would be singular at
#' this null (the effects of the shape parameters vanish when the moderation is 0). Curves that
#' change with p without such a parameter (e.g. a SHASH latent target correlated with the
#' predictor) are tested by their differences from the median (noted); curves constant by
#' construction are reported as such.
#' @param object A `birt_rtmb` fit.
#' @param p Quantile levels (at least two).
#' @param ... Unused.
#' @return A list of class `birt_qtest`: `test` (per curve: the tested `parameters`, `chisq`,
#'   `df`, `p.value`) and `differences` (per curve and level: `diff` = beta(p) - beta(median),
#'   `se`, `z`, `p.value`). With ELGM the statistic uses the posterior mean and covariance (a
#'   normal approximation).
#' @references Koenker, R., & Bassett, G. (1982). Robust tests for heteroscedasticity based on
#'   regression quantiles. *Econometrica, 50*(1), 43-61. \doi{10.2307/1912528}
#' @export
quantile_test <- function(object, p = object$control$p, ...) {
  if (!inherits(object, "birt_rtmb")) stop("quantile_test() needs a fit of birt(engine = \"rtmb\")", call. = FALSE)
  if (length(p) < 2) stop("give at least two quantile levels", call. = FALSE)
  q <- quantiles(object, p); V <- attr(q, "vcov")
  ref <- p[which.min(abs(p - 0.5))]
  nf <- natural_fun(object$S); nm <- names(nf(object$rtmb$par))
  wald <- function(sel) {                                     # Wald test of the selected natural parameters
    d <- delta(object, function(pp) nf(pp)[sel], cov = TRUE)
    if (!all(is.finite(d$V))) return(c(chisq = NA_real_, df = NA_real_))
    ev <- eigen(d$V, symmetric = TRUE); keep <- ev$values > 1e-10 * max(c(ev$values, 1e-300))
    c(chisq = sum((crossprod(ev$vectors[, keep, drop = FALSE], d$est))^2 / ev$values[keep]), df = sum(keep))
  }
  tests <- diffs <- list()
  for (g in split(seq_len(nrow(q)), paste(q$target, q$predictor))) {
    tg <- q$target[g[1]]; k <- q$predictor[g[1]]
    r <- g[q$p[g] == ref]; o <- setdiff(g, r)
    C <- matrix(0, length(o), nrow(q)); C[cbind(seq_along(o), o)] <- 1; C[, r] <- -1
    dv <- as.vector(C %*% q$est); W <- C %*% V %*% t(C); se <- sqrt(pmax(diag(W), 0))
    varies <- isTRUE(any(se > 1e-8 * max(1, abs(q$est[g]), na.rm = TRUE), na.rm = TRUE))
    sel <- intersect(c(sprintf("V(%s)|~|%s", tg, k), sprintf("%s|skew~|%s", tg, k)), nm)
    if (length(sel)) { w <- wald(sel); note <- "" }
    else if (varies) {
      ev <- eigen(W, symmetric = TRUE); keep <- ev$values > 1e-10 * max(ev$values)
      w <- c(chisq = sum((crossprod(ev$vectors[, keep, drop = FALSE], dv))^2 / ev$values[keep]), df = sum(keep))
      note <- "differences across p (no moderation parameter)"
    } else w <- c(chisq = NA_real_, df = 0)
    tests[[length(tests) + 1]] <- data.frame(target = tg, predictor = k,
      parameters = if (length(sel)) paste(sub("\\|~\\|", " ~ ", sub("\\|skew~\\|", " skew ~ ", sel)), collapse = ", ") else "",
      chisq = unname(w["chisq"]), df = as.integer(w["df"]), p.value = if (isTRUE(w["df"] > 0)) stats::pchisq(w["chisq"], w["df"], lower.tail = FALSE) else NA_real_,
      note = if (is.na(w["df"])) "no standard errors (see convergence(fit))" else if (w["df"] > 0) note else "constant in p by the model",
      stringsAsFactors = FALSE, row.names = NULL)
    if (!varies) next                                                  # constant curves: no differences to show
    z <- ifelse(se > 0, dv / se, NA_real_)
    diffs[[length(diffs) + 1]] <- data.frame(target = tg, predictor = k, p = q$p[o], diff = dv, se = ifelse(se > 0, se, NA_real_),
                                             z = z, p.value = 2 * stats::pnorm(-abs(z)), stringsAsFactors = FALSE)
  }
  structure(list(test = do.call(rbind, tests), differences = if (length(diffs)) do.call(rbind, diffs), ref = ref,
                 elgm = !is.null(object$rtmb$elgm)), class = "birt_qtest")
}

#' @rdname quantile_test
#' @param x A `birt_qtest` object.
#' @export
print.birt_qtest <- function(x, ...) {
  cat(sprintf("Change of quantile effects with p (Wald test of the scale / skew moderation by the predictor%s)\n", if (x$elgm) "; ELGM posterior, normal approximation" else ""))
  print_table(x$test, 3)
  if (!is.null(x$differences)) { cat(sprintf("\nDifferences from p = %s\n", format(x$ref))); print_table(x$differences, 3) }
  invisible(x)
}
