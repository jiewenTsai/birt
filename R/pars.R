# Parameter draws in model order, as matrices (draws x indicators). Used by the marginal
# likelihood, posterior predictive checks and item tables.

prep_pars <- function(fit, n_draws = NULL) {
  M <- as.matrix(fit$draws)
  nc <- coda::nchain(fit$draws); ni <- coda::niter(fit$draws)
  chain <- rep(seq_len(nc), each = ni)
  if (!is.null(n_draws) && n_draws < nrow(M)) {             # same number of draws from every chain
    per <- max(1, floor(n_draws / nc))
    rows <- unlist(lapply(seq_len(nc), function(c) (c - 1) * ni + unique(round(seq(1, ni, length.out = per)))))
    M <- M[rows, , drop = FALSE]; chain <- chain[rows]
  }
  sp <- fit$spec; pt <- fit$partab; S <- nrow(M); N <- fit$data$N
  val <- function(i) if (!is.na(pt$node[i])) M[, pt$node[i]] else rep(pt$fixed[i], S)
  lam <- function(f, vars) {
    out <- matrix(0, S, length(vars), dimnames = list(NULL, vars))
    if (is.na(f)) return(out)
    for (v in vars) { i <- which(pt$kind == "loading" & pt$lhs == f & pt$rhs == v); if (length(i)) out[, v] <- val(i) }
    out
  }
  cols <- function(vars, pat) if (!length(vars)) matrix(0, S, 0) else
    matrix(vapply(vars, function(v) M[, sprintf(pat, safe(v))], numeric(S)), nrow = S, dimnames = list(NULL, vars))
  f1 <- sp$factors[1]; f2 <- if (length(sp$factors) > 1) sp$factors[2] else NA_character_
  finfo <- function(f) {
    if (is.na(f)) return(NULL)
    fam <- sp$family$factor[[f]]
    i <- which(pt$lhs == f & pt$op == "~~" & pt$rhs == f)
    sc <- val(i)
    if (fam$family == "ald" && is.na(pt$node[i])) sc <- sc / sqrt(ald_var(fam$p))
    rr <- pt[pt$lhs == f & pt$kind == "beta", ]
    other <- setdiff(sp$factors, f)
    lin <- rr[!(rr$rhs %in% sp$factors), ]
    B <- matrix(if (nrow(lin)) M[, lin$node] else 0, S, nrow(lin), dimnames = list(NULL, lin$rhs))
    X <- if (nrow(lin)) matrix(unlist(fit$data[lin$rhs]), N, nrow(lin), dimnames = list(NULL, lin$rhs)) else matrix(0, N, 0)
    g <- if (length(other) && any(rr$rhs == other)) M[, rr$node[rr$rhs == other]] else rep(0, S)
    list(name = f, fam = fam, sc = sc, B = B, X = X, g = g)
  }
  rn <- which(pt$kind == "corr")
  # moderation: moderator values (N x Mz) and per draw the coefficients (S x Mz) of each loading,
  # shift, latent log SD and residual log SD (zero where there is no effect)
  zn <- sp$modvars %||% character()
  Z <- if (length(zn)) matrix(unlist(fit$data[zn]), N, length(zn), dimnames = list(NULL, zn)) else matrix(0, N, 0)
  coefs <- function(kind, lhs, rhs_fun) {
    out <- matrix(0, S, length(zn), dimnames = list(NULL, zn))
    for (z in zn) { i <- which(pt$kind == kind & pt$lhs == lhs & pt$rhs == rhs_fun(z)); if (length(i)) out[, z] <- val(i) }
    out
  }
  amod <- list(); bmod <- list(); smod <- list()
  if (length(zn)) {
    for (i in which(pt$kind == "loading")) amod[[paste(pt$lhs[i], pt$rhs[i])]] <- coefs("alpha", sprintf("E(%s)", pt$rhs[i]), function(z) paste0(z, ":", pt$lhs[i]))
    for (v in c(sp$binary, sp$ordinal, sp$cont)) bmod[[v]] <- coefs("delta", sprintf("E(%s)", v), identity)
    for (f in sp$factors) smod[[f]] <- coefs("psi", sprintf("V(%s)", f), identity)
    for (v in sp$cont) smod[[v]] <- coefs("kappa", sprintf("V(%s)", v), identity)
  }
  O <- sp$ordinal %||% character()
  thr <- lapply(stats::setNames(O, O), function(v) { K <- sp$ord$K[match(v, sp$ord$items)]
    matrix(vapply(seq_len(K - 1), function(k) M[, sprintf("b_%s[%d]", safe(v), k)], numeric(S)), nrow = S) })
  rat <- lapply(stats::setNames(O, O), function(v) { K <- sp$ord$K[match(v, sp$ord$items)]       # step ratios (tppcm)
    matrix(vapply(seq_len(K - 1), function(k) if (k > 1 && sp$ord$type[[v]] == "tppcm") M[, sprintf("r_%s[%d]", safe(v), k)] else rep(1, S), numeric(S)), nrow = S) })
  list(M = M, S = S, N = N, f1 = f1, f2 = f2, chain = chain, Z = Z, amod = amod, bmod = bmod, smod = smod,
       o1 = lam(f1, O), o2 = lam(f2, O), thr = thr, rat = rat, otype = sp$ord$type[O], Yo = if (length(O)) sp$ord$Y[, O, drop = FALSE] else NULL,
       a1 = lam(f1, sp$binary), a2 = lam(f2, sp$binary),
       l1 = lam(f1, sp$cont), l2 = lam(f2, sp$cont),
       d = cols(sp$binary, "d_%s"), xi = cols(sp$cont, "xi_%s"), sig = cols(sp$cont, "sigma_%s"),
       pc = vapply(sp$cont, function(v) { f <- sp$family$ind[[v]]; if (is_q(f)) f$p else NA_real_ }, 0),
       ct = vapply(sp$cont, function(v) sp$family$ind[[v]]$family, ""),
       F1 = finfo(f1), F2 = finfo(f2),
       r = if (length(rn)) M[, pt$node[rn]] else rep(0, S))
}

# moderation at draw s: exp(alpha'z) of a loading, delta'z of a shift, exp(psi'z) of an SD (N-vectors;
# 1 / 0 / 1 without moderation)
mod_mult <- function(P, s, key) { a <- P$amod[[key]]; if (is.null(a) || !any(a[s, ] != 0)) 1 else exp(drop(P$Z %*% a[s, ])) }
mod_shift <- function(P, s, v) { b <- P$bmod[[v]]; if (is.null(b) || !any(b[s, ] != 0)) 0 else drop(P$Z %*% b[s, ]) }
mod_sd <- function(P, s, x) { a <- P$smod[[x]]; if (is.null(a) || !any(a[s, ] != 0)) 1 else exp(drop(P$Z %*% a[s, ])) }

# conditional mean of a factor from covariates, N-vector
mean_cov <- function(Fi, s, N) if (ncol(Fi$B)) drop(Fi$X %*% Fi$B[s, ]) else rep(0, N)
# contribution of the latent predictor th (N-vector) to the mean of the second factor
mean_path <- function(Fi, s, th) Fi$g[s] * th
