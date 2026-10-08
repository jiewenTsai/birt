# Approximate Bayesian inference with the RTMB engine (rtmb_control(method = "elgm")) ---------
#
# The model is an extended latent Gaussian model (Stringer, Brown & Stafford, 2023) with the
# person latents z as the latent Gaussian field W and every other parameter theta in the outer
# layer:
#   p(theta | y) ~ p(y | theta) p(theta),  p(y | theta) = prod_i int p(y_i | z_i, theta) phi(z_i) dz_i
# The inner integrals are per-person AGHQ (or Laplace). The outer integral uses a mixed AGHQ
# grid at the posterior mode: k_hyper nodes on the `hyper` directions, which come first in the
# Cholesky factor of the inverse Hessian so that their node spread is marginal; the other
# directions sit at their Gaussian conditional modes and add their conditional variance
# (pilot_20260928: R/elgm.R, findings lesson 4). Item parameters are not put into W (pilot
# lesson 1: the inner mode would be a penalized JML with incidental-parameter bias).

# ---- priors on the internal (working) parameters -------------------------------------------
# Each term: working parameter `par[idx]`, transform to the natural scale ("id", "exp", "tanh"),
# the natural-scale log density, and the name of the natural quantity in the estimates table.
rtmb_priors <- function(S, dp) {
  terms <- list(); txt <- character()
  add <- function(par, idx, tr, prior, key, show = prior) {
    dens <- prior_density(prior)
    lo <- switch(tr, exp = 0, tanh = -1, -Inf); hi <- switch(tr, tanh = 1, Inf)
    if (dens$lo < lo || dens$hi > hi) {            # truncate to the parameter's range (e.g. dnorm on a positive loading)
      if (grepl("T\\(", prior)) stop(sprintf("the prior %s of %s allows values outside (%g, %g)", prior, key, lo, hi), call. = FALSE)
      prior <- sprintf("%s T(%s,%s)", prior, if (is.finite(lo)) lo else "", if (is.finite(hi)) hi else "")
      dens <- prior_density(prior)
    } else if (dens$lo > lo && tr == "id") stop(sprintf("the prior %s of %s excludes negative values; this parameter is unrestricted", prior, key), call. = FALSE)
    terms[[length(terms) + 1]] <<- list(par = par, idx = idx, tr = tr, f = dens$f)
    txt[key] <<- show
  }
  ld <- S$ld; seen <- list(); hier <- list()
  for (r in seq_len(nrow(ld))) {
    if (ld$kind[r] == "fixed") next
    key <- sprintf("%s|=~|%s", ld$lhs[r], ld$rhs[r])
    pr <- if (ld$prior[r] != "") ld$prior[r] else if (ld$first[r]) dp$loading else dp$cross
    if (pr == "ssp") stop("prior(\"ssp\") needs engine = \"jags\"; with method = \"elgm\" use a small-variance normal prior, e.g. prior(\"dnorm(0, 100)\")", call. = FALSE)
    par <- if (ld$kind[r] == "pos") "lpos" else "lreal"
    # loadings sharing a label share one parameter (and the prior of its first line)
    id <- paste(par, ld$idx[r])
    if (id %in% names(seen)) { txt[key] <- seen[[id]]; next }
    seen[[id]] <- pr
    if (pr == "hier") {                                   # joint hierarchical prior of the positive loadings of a latent variable
      if (par != "lpos") stop("the hierarchical loading prior (\"hier\") is for positive loadings", call. = FALSE)
      hier[[ld$lhs[r]]] <- c(hier[[ld$lhs[r]]], ld$idx[r]); seen[[id]] <- txt[key] <- sprintf("hier(%s)", ld$lhs[r]); next
    }
    add(par, ld$idx[r], if (par == "lpos") "exp" else "id", pr, key)
  }
  for (i in seq_len(S$nF)) {
    fi <- S$f[[i]]; f <- fi$name
    if (i == 2 && fi$path) add("gam", 1, "id", if (fi$path_prior != "") fi$path_prior else dp$beta, sprintf("%s|~|%s", f, S$fs[1]))
    for (k in seq_along(fi$covs)) {
      pr <- if (fi$cov_prior[k] != "") fi$cov_prior[k] else dp$beta
      if (pr == "ssp") stop("prior(\"ssp\") needs engine = \"jags\"", call. = FALSE)
      add(paste0("beta", i), k, "id", pr, sprintf("%s|~|%s", f, fi$covs[k]))
    }
    plain <- !fi$shash && !length(fi$scale)
    if (fi$free_sd) add(paste0("lsd", i), 1, "exp", dp$factor_sd, if (plain) sprintf("%s|~~|%s", f, f) else sprintf("%s|sd|", f),
                        show = if (plain) paste("sd", dp$factor_sd) else dp$factor_sd)
    sprior <- function(pr, z) { x <- unname(pr[z]); if (length(x) && !is.na(x) && nzchar(x)) x else dp$moderation }
    for (k in seq_along(fi$scale)) add(paste0("alpha", i), k, "id", sprior(fi$scale_prior, fi$scale[k]), sprintf("V(%s)|~|%s", f, fi$scale[k]))
    if (fi$shash) {
      add("eps2", 1, "id", dp$skew, sprintf("%s|skew|", f))
      for (k in seq_along(fi$skew)) add("eta2", k, "id", dp$moderation, sprintf("%s|skew~|%s", f, fi$skew[k]))
      add("ldl2", 1, "exp", dp$tail, sprintf("%s|tail|", f))
    }
  }
  if (S$cor) add("atr", 1, "tanh", dp$cor, sprintf("%s|~~|%s", S$fs[1], S$fs[2]))
  for (j in seq_along(S$B)) add("d", j, "id", dp$intercept_binary, sprintf("%s|~1|", S$B[j]))
  o <- 0
  for (j in seq_along(S$O)) {                                            # thresholds: first, then the positive steps
    ix <- S$tix[[j]]
    if (S$otype[j] != "grm") {
      for (k in seq_along(ix$pb)) add("pb", ix$pb[k], "id", dp$threshold, sprintf("%s|thr|t%d", S$O[j], k))
      for (k in seq_along(ix$lr)) add("lr", ix$lr[k], "exp", dp$step_ratio, sprintf("%s|stepa|a%d", S$O[j], k + 1), show = sprintf("loading x %s", dp$step_ratio))
      next
    }
    add("t1", ix$t1, "id", dp$threshold, sprintf("%s|thr|t1", S$O[j]))
    for (k in seq_along(ix$linc)) add("linc", ix$linc[k], "exp", dp$threshold_step, sprintf("%s|thr|t%d", S$O[j], k + 1))
  }
  if (!is.null(S$mod)) for (pn in c("af", "ac", "bf", "bc")) for (i in unique(S$mod$idx[S$mod$par == pn])) {
    m <- S$mod[S$mod$par == pn & S$mod$idx == i, ][1, ]
    pr <- if (nzchar(m$prior)) m$prior else if (m$kind == "alpha") dp$moderation else dp$beta
    for (r in which(S$mod$par == pn & S$mod$idx == i)) {
      mm <- S$mod[r, ]; key <- sprintf("E(%s)|~|%s", mm$target, if (mm$kind == "alpha") paste0(mm$mod, ":", mm$factor) else mm$mod)
      if (r == which(S$mod$par == pn & S$mod$idx == i)[1]) add(pn, i, "id", pr, key) else txt[key] <- pr
    }
  }
  for (j in seq_along(S$C)) {
    v <- S$C[j]; kr <- which(S$kap$j == j); s <- match(j, S$sh)
    add("xi", j, "id", dp$intercept, sprintf("%s|~1|", v))
    plain <- is.na(s) && !length(kr)
    add("lsig", j, "exp", dp$resid_sd, if (plain) sprintf("%s|~~|%s", v, v) else sprintf("%s|sd|", v),
        show = if (plain) paste("sd", dp$resid_sd) else dp$resid_sd)
    for (r in kr) { sp0 <- S$famc[[j]]$scale_prior %||% character(); x <- unname(sp0[S$kap$m[r]])
      add("kap", r, "id", if (length(x) && !is.na(x) && nzchar(x)) x else dp$moderation, sprintf("V(%s)|~|%s", v, S$kap$m[r])) }
    if (!is.na(s)) {
      add("eps", s, "id", dp$skew, sprintf("%s|skew|", v))
      for (r in which(S$eta$j == j)) add("eta", r, "id", dp$moderation, sprintf("%s|skew~|%s", v, S$eta$m[r]))
      add("ldl", s, "exp", dp$tail, sprintf("%s|tail|", v))
    }
  }
  # log prior of the working parameters: natural-scale density + log Jacobian
  lp <- function(p) {
    out <- 0
    for (t in terms) {
      w <- p[[t$par]][t$idx]
      out <- out + switch(t$tr, id = t$f(w), exp = t$f(exp(w)) + w, tanh = { r <- tanh(w); t$f(r) + log(1 - r * r) })
    }
    for (h in hier) out <- out + hier_lp(p$lpos[h])      # density of log lambda: no Jacobian
    out
  }
  list(lp = lp, txt = txt, n = length(terms) + length(unlist(hier)))
}

# numeric vector (names as in obj$par) -> parameter list with the skeleton's structure
vec2list <- function(x, skel) {
  sp <- split(unname(x), factor(names(x), levels = unique(names(x))))
  for (n in names(skel)) skel[[n]] <- if (n %in% names(sp)) sp[[n]] else skel[[n]][0]
  skel
}

# per-person AGHQ log-likelihood (numeric, without the prior) and the EAP moments
aghq_value <- function(S, p, cn, k) {
  ND <- node_data(S, cn, k); N <- S$N
  Fl <- lat_values(p, S, ND$z1, ND$z2, ND$rows)
  L <- matrix(ind_ll(p, S, Fl, ND$rows) + ND$lw, N)
  ls <- lse(L)
  W <- exp(L - ls)
  mom <- lapply(seq_len(S$nF), function(i) { Fm <- matrix(Fl[[i]], N); cbind(m = rowSums(W * Fm), m2 = rowSums(W * Fm^2)) })
  list(ll = sum(ls + ND$cst), mom = mom)
}

hyper_names <- function(hyper) {
  key <- list(cor = "atr", path = "gam", sd = c("lsd1", "lsd2"), shape = c("eps2", "ldl2", "eta2", "alpha1", "alpha2"),
              cross = "lreal", resid = "lsig")
  unique(unlist(lapply(hyper, function(h) key[[h]] %||% h)))
}

# outer mixed-grid AGHQ at the posterior mode r (from fit_rtmb with S$prior set)
elgm_outer <- function(S, r, ctrl, progress = TRUE) {
  say <- function(...) if (progress) message(...)
  t0 <- Sys.time()
  x <- r$x; nx <- length(x); skel <- r$par
  free <- !r$ab
  V <- r$vcov
  ok <- free & is.finite(diag(V)) & diag(V) > 0
  if (any(free & !ok)) warning("ELGM: the Hessian at the posterior mode is not positive definite in some directions; they are fixed at the mode", call. = FALSE)
  hyper <- ctrl$hyper
  if (identical(hyper, "auto")) hyper <- if (any(names(x) %in% hyper_names("shape"))) "shape" else c("cor", "path", "sd")
  hn <- hyper_names(hyper)
  bad <- setdiff(hn, c(names(x), hyper_names(c("cor", "path", "sd", "shape", "cross", "resid"))))
  if (length(bad)) stop("rtmb_control(hyper = ): unknown parameters: ", paste(bad, collapse = ", "), "; see names(fit$rtmb$x)", call. = FALSE)
  q <- which(names(x) %in% hn & ok)
  if (identical(ctrl$hyper, "auto") && length(q) > 3) {
    # at most 3 quadrature directions under "auto" (k_hyper = 5: 125 instead of 625+ nodes): the
    # SHASH shape first, then the moderation effects, then correlation / path / SDs; the rest are
    # Gaussian directions (results/birt_groupC/84_elgm_hyper_dims.txt: 3 of 4 directions changed the
    # posterior means by <= 0.034 SD, the SDs by <= 2.1% and log p(y) by 0.02, at 1/5 of the outer time)
    pri <- c("eps2", "ldl2", "eta2", "alpha2", "alpha1", "atr", "gam", "lsd2", "lsd1")
    keep <- sort(q[order(match(names(x)[q], pri), q)][1:3]); drop <- setdiff(q, keep)
    message(sprintf("ELGM: hyper = \"auto\" uses outer quadrature on at most 3 hyperparameters (%s); %s %s treated as Gaussian. For quadrature on all %d give rtmb_control(hyper = c(%s)) (%d^%d nodes)",
                    paste(names(x)[keep], collapse = ", "), paste(names(x)[drop], collapse = ", "), if (length(drop) == 1) "is" else "are",
                    length(q), paste(sprintf("\"%s\"", unique(names(x)[q])), collapse = ", "), c(1, 11, 9, 5, 5, 3, 3)[min(length(q), 6) + 1], length(q)))
    q <- keep
  }
  nq <- length(q)
  k <- ctrl$k_hyper %||% c(1, 11, 9, 5, 5, 3, 3)[min(nq, 6) + 1]
  if (nq && k^nq > 2000) stop(sprintf("ELGM: %d^%d = %d outer nodes; choose fewer hyperparameters (hyper = ) or a smaller k_hyper", k, nq, k^nq), call. = FALSE)
  rest <- setdiff(which(ok), q); ordf <- c(q, rest)
  Vo <- V[ordf, ordf, drop = FALSE]; Vo <- (Vo + t(Vo)) / 2
  Lc <- t(chol(Vo))
  gh <- statmod::gauss.quad.prob(k, "normal")
  Zq <- if (nq) as.matrix(expand.grid(rep(list(gh$nodes), nq))) else matrix(0, 1, 0)
  lwq <- if (nq) rowSums(log(as.matrix(expand.grid(rep(list(gh$weights), nq))))) - rowSums(stats::dnorm(Zq, log = TRUE)) else 0
  G <- nrow(Zq)
  # conditional covariance of the Gaussian directions given the quadrature ones
  Sc <- matrix(0, nx, nx)
  Sc[ok, ok] <- V[ok, ok] - if (nq) V[ok, q, drop = FALSE] %*% solve(V[q, q, drop = FALSE], V[q, ok, drop = FALSE]) else 0
  # one Laplace object for the inner modes (and the inner Laplace value)
  z0 <- list(z1 = rep(0, S$N)); if (S$nF == 2) z0$z2 <- rep(0, S$N)
  objL <- RTMB::MakeADFun(make_joint(S), c(skel, z0), random = re_names(S), silent = TRUE)
  kin <- r$k
  # conditional modes of the Gaussian directions at each node: optimized on the AGHQ (or
  # Laplace) objective with the quadrature centres of the mode (taped once)
  if (ctrl$cond_modes && nq && length(rest)) {
    objL$fn(x)
    objC <- if (r$method == "aghq") RTMB::MakeADFun(make_aghq(S, node_data(S, centers_at(S, objL), kin)), skel, silent = TRUE) else objL
    lo <- rtmb_bounds(x, "lo")[rest]; up <- rtmb_bounds(x, "up")[rest]
    # Newton steps with the conditional covariance of the mode (H_rr^-1 = Sc[rest, rest]) and step halving
    Srr <- Sc[rest, rest, drop = FALSE]
    cmode <- function(xg) {
      xg[rest] <- pmin(pmax(xg[rest], lo), up)
      fv <- tryCatch(objC$fn(xg), error = function(e) NA)
      if (!is.finite(fv)) return(xg)
      for (it in seq_len(8)) {
        g <- tryCatch(objC$gr(xg)[rest], error = function(e) NULL)
        if (is.null(g) || any(!is.finite(g))) break
        step <- -as.vector(Srr %*% g); dec <- -sum(g * step)
        if (dec < 1e-6) break
        acc <- FALSE
        for (t in c(1, 0.5, 0.25, 0.125)) {
          y <- xg; y[rest] <- pmin(pmax(xg[rest] + t * step, lo), up)
          fy <- tryCatch(objC$fn(y), error = function(e) NA)
          if (is.finite(fy) && fy < fv) { xg <- y; acc <- fv - fy > 1e-7; fv <- fy; break }
        }
        if (!acc) break
      }
      xg
    }
  } else cmode <- identity
  say(sprintf("ELGM: outer AGHQ, %d nodes (%s)", G,
              if (nq) sprintf("k = %d on %s; Gaussian in %d other directions", k, paste(unique(names(x)[q]), collapse = ", "), length(rest))
              else "posterior mode + Laplace"))
  # node at standardized position z (quadrature directions), Gaussian directions at their conditional modes
  node_x <- function(z) {
    xg <- x
    if (nq) xg[ordf] <- xg[ordf] + as.vector(Lc %*% c(z, numeric(length(rest))))
    if (any(z != 0)) cmode(xg) else xg
  }
  eval_x <- function(xg) {
    v <- tryCatch(objL$fn(xg), error = function(e) NA)
    if (!is.finite(v)) return(list(lpost = -Inf, ll = NA, mom = NULL))
    p <- vec2list(xg, skel)
    A <- aghq_value(S, p, centers_at(S, objL), kin)
    list(lpost = if (r$method == "laplace") -v else A$ll + S$prior$lp(p), ll = A$ll, mom = A$mom)
  }
  # asymmetric scaling of each quadrature direction (as INLA's grid): from the drop of the log
  # posterior at z = +-2, d = sqrt(2 / drop) (1 for a Gaussian posterior), bounded to [1/3, 3]
  e0 <- eval_x(x)
  dsc <- matrix(1, max(nq, 1), 2, dimnames = list(NULL, c("neg", "pos")))
  if (nq && ctrl$cond_modes) for (i in seq_len(nq)) for (sg in 1:2) {
    z <- numeric(nq); z[i] <- c(-2, 2)[sg]
    drop <- e0$lpost - eval_x(node_x(z))$lpost
    dsc[i, sg] <- if (is.finite(drop)) min(max(sqrt(2 / max(drop, 1e-3)), 1 / 3), 3) else 1 / 3
  }
  Zs <- Zq
  if (nq) for (i in seq_len(nq)) Zs[, i] <- Zq[, i] * ifelse(Zq[, i] < 0, dsc[i, 1], dsc[i, 2])
  ljac <- if (nq) rowSums(vapply(seq_len(nq), function(i) log(ifelse(Zq[, i] < 0, dsc[i, 1], dsc[i, 2])), numeric(G))) else 0
  logw <- lwq + ljac + sum(log(diag(Lc))) + length(rest) * 0.5 * log(2 * pi)
  X <- matrix(x, G, nx, byrow = TRUE, dimnames = list(NULL, names(x)))
  lpost <- numeric(G); ll <- numeric(G); mom <- vector("list", G)
  for (g in seq_len(G)) {
    if (nq && any(Zq[g, ] != 0)) { X[g, ] <- node_x(Zs[g, ]); ev <- eval_x(X[g, ]) } else ev <- e0
    lpost[g] <- ev$lpost; ll[g] <- ev$ll; mom[g] <- list(ev$mom)
  }
  if (all(!is.finite(lpost))) stop("ELGM: the log posterior could not be evaluated at the outer nodes", call. = FALSE)
  lw <- logw + lpost; mx <- max(lw); lnc <- mx + log(sum(exp(lw - mx)))
  prob <- exp(lw - lnc)
  centre <- which.max(-rowSums(Zq^2))
  excess <- max(lpost) - lpost[centre]
  if (nq && excess > 0.5)
    warning(sprintf("ELGM: an outer node has a log posterior %.2f above the mode; the posterior may be multimodal or the mode inaccurate", excess), call. = FALSE)
  # grid coverage: mass on the two outermost nodes of each quadrature direction against its
  # value for a Gaussian posterior (2 w_min of the k-node rule)
  ec <- elgm_edge(prob, Zq, gh)
  edge <- if (nq) max(ec$mass) else 0; edge0 <- if (nq) ec$gauss[1] else 0
  if (nq && any(ec$wide)) {
    i <- which.max(ec$mass / ec$limit)
    warning(sprintf("ELGM: the posterior is too wide for the outer grid along %s: %.1f%% of the mass is on its two outermost nodes (Gaussian posterior: %.1f%%; limit: %.1f%%; k_hyper = %d), so the grid may cut off a tail of the posterior: check whether the parameter is weakly determined (e.g. a variance near 0) and refit with a larger k_hyper",
                    names(x)[q][i], 100 * ec$mass[i], 100 * ec$gauss[i], 100 * ec$limit[i], k), call. = FALSE)
  }
  # EAP scores: mixtures over the outer nodes
  okg <- which(prob > 0 & !vapply(mom, is.null, TRUE))
  sc <- lapply(seq_len(S$nF), function(i) {
    m <- Reduce(`+`, lapply(okg, function(g) prob[g] * mom[[g]][[i]][, "m"])) / sum(prob[okg])
    m2 <- Reduce(`+`, lapply(okg, function(g) prob[g] * mom[[g]][[i]][, "m2"])) / sum(prob[okg])
    list(mean = m, sd = sqrt(pmax(m2 - m^2, 0)))
  })
  list(nodes = X, prob = prob, logpost = lpost, loglik = ll, lognc = lnc, hyper = unique(names(x)[q]), nq = nq, k = k,
       n_gauss = length(rest), Sc = Sc, scaling = if (nq) cbind(par = names(x)[q], as.data.frame(dsc)), scores = sc, ess = 1 / sum(prob^2), excess = excess, edge = edge, edge0 = edge0, coverage = if (nq) cbind(par = names(x)[q], ec),
       cond_modes = ctrl$cond_modes && nq > 0 && length(rest) > 0,
       inner = r$method, k_inner = kin, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# Coverage of the outer grid (in the spirit of the grid exploration of INLA, Rue, Martino &
# Chopin, 2009, sec. 6.5: the grid has to cover the region where the posterior has mass).
# Along quadrature direction i, m_i is the posterior mass on the two outermost nodes
# (|z| = z_max). For a Gaussian posterior s times wider than the rule assumes (s = 1 after a
# correct scaling), the k-node rule puts
#   m(s; k) = sum_{|z_n| = z_max} w_n r_s(z_n) / sum_n w_n r_s(z_n),  r_s(z) = phi(z / s) / (s phi(z))
# there (s = 1: 2 w_min, i.e. 1/3 for k = 3, .023 for k = 5, 2.4e-7 for k = 11). A direction
# is flagged when m_i > m(1.5; k), i.e. its outermost nodes carry as much mass as under a
# Gaussian posterior 1.5 times wider than the grid assumes (the same width for every k), and
# m_i > .025 (for k >= 9 the outermost nodes are beyond 4.5 SD, where heavier than Gaussian
# tails of a well covered posterior already exceed m(1.5; k) = .002-.009). Both thresholds are
# heuristic choices (not from the literature); they flag a Gaussian posterior wider than
# 1.5 (k = 3, 5), 1.76 (k = 9) or 2.09 (k = 11) times the assumed width.
gauss_edge_mass <- function(gh, s) {
  w <- gh$weights * stats::dnorm(gh$nodes / s) / (s * stats::dnorm(gh$nodes))
  sum(w[abs(gh$nodes) >= max(abs(gh$nodes)) - 1e-8]) / sum(w)
}
elgm_edge <- function(prob, Zq, gh, wide = 1.5) {
  if (!ncol(Zq)) return(data.frame(mass = numeric(), gauss = numeric(), limit = numeric(), wide = logical()))
  zmax <- max(abs(gh$nodes)); lim <- gauss_edge_mass(gh, wide)
  m <- vapply(seq_len(ncol(Zq)), function(i) sum(prob[abs(Zq[, i]) >= zmax - 1e-8]), 0)
  data.frame(mass = m, gauss = gauss_edge_mass(gh, 1), limit = pmax(lim, 0.025), wide = m > pmax(lim, 0.025))
}

# posterior summaries of g(theta) under the outer quadrature (mean, SD and a 95% interval from
# the quadrature moments: Cornish-Fisher with the skewness on the scale tr = "log" / "atanh" / "id")
elgm_summary <- function(E, G, J, tr) {
  pr <- E$prob
  vc <- pmax(rowSums((J %*% E$Sc) * J), 0)
  m <- as.vector(G %*% pr)
  sdv <- sqrt(pmax(as.vector((G - m)^2 %*% pr) + vc, 0))
  h <- function(v, t) switch(t, log = log(pmax(v, 1e-300)), atanh = atanh(pmax(pmin(v, 1 - 1e-12), -1 + 1e-12)), v)
  hi <- function(v, t) switch(t, log = exp(v), atanh = tanh(v), v)
  dh <- function(v, t) switch(t, log = 1 / v, atanh = 1 / (1 - v^2), 1)
  lo <- up <- numeric(length(m))
  for (i in seq_along(m)) {
    t <- tr[i]
    if (t != "id" && (t == "log" && any(G[i, ] <= 0) || t == "atanh" && any(abs(G[i, ]) >= 1))) t <- "id"
    H <- h(G[i, ], t); mh <- sum(pr * H)
    vh <- sum(pr * (H - mh)^2) + dh(m[i], t)^2 * vc[i]
    g1 <- if (vh > 0) max(min(sum(pr * (H - mh)^3) / vh^1.5, 1), -1) else 0
    cf <- function(z) hi(mh + sqrt(vh) * (z + g1 * (z^2 - 1) / 6), t)
    lo[i] <- cf(-1.959964); up[i] <- cf(1.959964)
  }
  list(est = m, se = sdv, lo = lo, hi = up)
}
