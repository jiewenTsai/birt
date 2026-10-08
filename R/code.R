# JAGS code generation ---------------------------------------------------------------
#
# Returns the model code, the parameter table (one row per model parameter, with the
# JAGS node that holds it) and the nodes to monitor.
#
# Node names: lam_<f>_<v> loadings (lab_<label> for labelled ones), d_<v> / xi_<v>
# (labb_<label> for labelled regression coefficients)
# intercepts, sigma_<v> residual SD (ALD: scale), sd_<f> / scale_<f> factor SD / ALD
# scale, r_<f1>_<f2> factor correlation, beta_<f>_<x> regressions, incl_* spike-and-slab indicators,
# e_* ALD mixing variables (not monitored).

# Lines for a variable X (data name xn, node suffix sx) with location mu and an ALD family.
# sc: scale expression (node name or number); sdn: the implied SD (node name, or a number when
# fixed); the scale is derived from it. Returns list(person, global).
q_lines <- function(xn, sx, mu, fam, sc, sdn) {
  k <- ald_k(fam$p)
  list(person = c(sprintf("    e_%s[j] ~ dexp(1 / %s)", sx, sc),
                  sprintf("    %s[j] ~ dnorm(%s%s, 1 / (%s * %s * e_%s[j]))", xn, mu,
                          if (k[1] != 0) sprintf(" + %s * e_%s[j]", num(k[1]), sx) else "",
                          num(k[2]), sc, sx)),
       global = sprintf("  %s <- %s / %s", sc, sdn, num(sqrt(ald_var(fam$p)))))
}

build_code <- function(spec) {
  ld <- spec$ld; fs <- spec$factors; B <- spec$binary; C <- spec$cont; O <- spec$ordinal; rg <- spec$reg
  fam <- spec$family; dp <- spec$dp %||% dpriors(); ssp_p <- spec$ssp_p
  mo <- spec$mod %||% mod_empty()
  lnode <- function(f, v) sprintf("lam_%s_%s", f, safe(v))
  bnode <- function(f, x) sprintf("beta_%s_%s", f, safe(x))
  rows <- list()
  add <- function(lhs, op, rhs, node = NA_character_, fixed = NA_real_, prior = "", kind = "", incl = NA_character_, group = NA_character_)
    rows[[length(rows) + 1]] <<- data.frame(lhs = lhs, op = op, rhs = rhs, node = node, fixed = fixed,
                                             prior = prior, kind = kind, incl = incl, group = group, stringsAsFactors = FALSE)
  # moderation effects: one node per free / spike-and-slab cell, one per label
  mnode <- function(r) if (r$type == "common") sprintf("%sc_%s", if (r$kind == "alpha") "alm" else "dlt", safe(r$label))
                       else sprintf("%s_%s_%s%s", if (r$kind == "alpha") "alm" else "dlt", safe(r$target),
                                    if (r$kind == "alpha") paste0(safe(r$factor), "_") else "", safe(r$mod))
  msum <- function(rr) paste(vapply(seq_len(nrow(rr)), function(i) sprintf("%s * %s[j]", mnode(rr[i, ]), rr$mod[i]), ""), collapse = " + ")
  act <- mo[mo$type != "none", ]
  # loading of v on f at person j: lam * exp(alpha'z) when moderated
  lexpr <- function(f, v) if (any(act$kind == "alpha" & act$target == v & act$factor == f)) sprintf("lm_%s_%s[j]", f, safe(v)) else lnode(f, v)
  term <- function(v) {
    r <- ld[ld$rhs == v & !ld$zero, ]
    if (!nrow(r)) return("0")
    paste(sprintf("%s * %s[j]", vapply(seq_len(nrow(r)), function(i) lexpr(r$lhs[i], v), ""), r$lhs), collapse = " + ")
  }
  has_shift <- function(v) any(act$kind == "beta" & act$target == v)
  smods <- function(f) setdiff(f$scale %||% character(), fs)                  # observed scale moderators (V() ~ z)
  sexpr <- function(base, f, pre, nm) if (!length(smods(f))) base else sprintf("%s_%s_j[j]", pre, nm)

  # ---- measurement part -------------------------------------------------------------
  lik <- character()
  for (r in seq_len(nrow(ld))) {                                                     # moderated loadings
    v <- ld$rhs[r]; f <- ld$lhs[r]; aa <- act[act$kind == "alpha" & act$target == v & act$factor == f, ]
    if (nrow(aa) && !ld$zero[r]) lik <- c(lik, sprintf("    lm_%s_%s[j] <- %s * exp(%s)", f, safe(v), lnode(f, v), msum(aa)))
  }
  for (v in c(B, O, C)) if (has_shift(v)) lik <- c(lik, sprintf("    dl_%s[j] <- %s", safe(v), msum(act[act$kind == "beta" & act$target == v, ])))
  shift <- function(v) {                                                             # on the latent scale: - loading * delta
    if (!has_shift(v)) return("")
    r <- ld[ld$rhs == v & !ld$zero, ]
    sprintf(" - %s * dl_%s[j]", lexpr(r$lhs[1], v), safe(v))
  }
  for (v in B) lik <- c(lik, sprintf("    logit(p_%s[j]) <- d_%s + %s%s", safe(v), safe(v), term(v), shift(v)),
                        sprintf("    %s[j] ~ dbern(p_%s[j])", v, safe(v)))
  for (i in seq_along(O)) {                                                         # P(y > k) = logit^-1{a (f - delta) - a b_k}
    v <- O[i]; sv <- safe(v); K <- spec$ord$K[match(v, spec$ord$items)]; r <- ld[ld$rhs == v & !ld$zero, ]
    a <- lexpr(r$lhs, v); pj <- grepl("[j]", a, fixed = TRUE); ty <- spec$ord$type[[v]]
    if (ty != "grm") {                                                               # adjacent-category logits: weights exp(cumulative sums)
      ak <- if (ty == "tppcm") sprintf("%s * r_%s[k - 1]", a, sv) else a
      lik <- c(lik, sprintf("    lp_%s[j, 1] <- 0", sv),
               sprintf("    for (k in 2:%d) { lp_%s[j, k] <- lp_%s[j, k - 1] + %s * (%s[j]%s - b_%s[k - 1]) }", K, sv, sv, ak, r$lhs,
                       if (has_shift(v)) sprintf(" - dl_%s[j]", sv) else "", sv),
               sprintf("    mx_%s[j] <- max(lp_%s[j, 1:%d])", sv, sv, K),
               sprintf("    for (k in 1:%d) { ep_%s[j, k] <- exp(lp_%s[j, k] - mx_%s[j]) }", K, sv, sv, sv),
               sprintf("    %s[j] ~ dcat(ep_%s[j, 1:%d])", v, sv, K))
      next
    }
    if (pj) lik <- c(lik, sprintf("    for (k in 1:%d) { cut_%s[j, k] <- %s * b_%s[k] }", K - 1, sv, a, sv))
    lik <- c(lik, sprintf("    %s[j] ~ dordered.logit(%s * (%s[j]%s), %s)", v, a, r$lhs,
                          if (has_shift(v)) sprintf(" - dl_%s[j]", sv) else "",
                          if (pj) sprintf("cut_%s[j, 1:%d]", sv, K - 1) else sprintf("cut_%s[1:%d]", sv, K - 1)))
  }
  qglob <- character()
  for (v in C) {
    f <- fam$ind[[v]]; sv <- safe(v)
    mu <- sprintf("xi_%s + %s%s", sv, term(v), if (has_shift(v)) sprintf(" + dl_%s[j]", sv) else "")
    if (length(smods(f))) {
      if (f$family != "normal") stop("V(", v, ") ~ z needs a normal residual with engine = \"jags\"", call. = FALSE)
      lik <- c(lik, sprintf("    sigma_%s_j[j] <- sigma_%s * exp(%s)", sv, sv, paste(sprintf("kap_%s_%s * %s[j]", sv, safe(smods(f)), smods(f)), collapse = " + ")))
    }
    if (f$family == "normal") lik <- c(lik, sprintf("    %s[j] ~ dnorm(%s, 1 / pow(%s, 2))", v, mu, sexpr(sprintf("sigma_%s", sv), f, "sigma", sv)))
    else {
      ql <- q_lines(v, sv, mu, f, sprintf("sigma_%s", sv), sprintf("sdr_%s", sv))
      lik <- c(lik, ql$person); qglob <- c(qglob, ql$global)
    }
  }

  # ---- structural part ----------------------------------------------------------------
  persons <- character()
  for (f in fs) {
    r <- rg[rg$lhs == f, ]
    terms <- character()
    for (i in seq_len(nrow(r))) {
      x <- r$x[i]
      terms <- c(terms, sprintf("%s * %s[j]", bnode(f, x), x))
    }
    mu <- if (length(terms)) {
      persons <- c(persons, sprintf("    mu_%s[j] <- %s", f, paste(terms, collapse = " + ")))
      sprintf("mu_%s[j]", f)
    } else "0"
    ff <- fam$factor[[f]]
    sdf <- if (is.na(spec$fsd[[f]])) sprintf("sd_%s", f) else num(spec$fsd[[f]])
    if (length(smods(ff))) {
      if (is_q(ff)) stop("V(", f, ") ~ z cannot be combined with ald() for ", f, call. = FALSE)
      persons <- c(persons, sprintf("    sd_%s_j[j] <- %s * exp(%s)", f, sdf, paste(sprintf("psi_%s_%s * %s[j]", f, safe(smods(ff)), smods(ff)), collapse = " + ")))
      sdf <- sprintf("sd_%s_j[j]", f)
    }
    if (is_q(ff)) {
      ql <- q_lines(f, f, mu, ff, sprintf("scale_%s", f), sdf)
      persons <- c(persons, ql$person); qglob <- c(qglob, ql$global)
    } else if (f == fs[2] && isTRUE(spec$cov_free)) {
      f1 <- fs[1]; sd1 <- if (length(smods(fam$factor[[f1]]))) sprintf("sd_%s_j[j]", f1) else if (is.na(spec$fsd[[f1]])) sprintf("sd_%s", f1) else num(spec$fsd[[f1]])
      mu1 <- if (any(rg$lhs == f1)) sprintf("mu_%s[j]", f1) else "0"
      persons <- c(persons, sprintf(
        "    %s[j] ~ dnorm(%s + r_%s_%s * %s / %s * (%s[j] - %s), 1 / (pow(%s, 2) * (1 - r_%s_%s * r_%s_%s)))",
        f, mu, f1, f, sdf, sd1, f1, mu1, sdf, f1, f, f1, f))
    } else persons <- c(persons, sprintf("    %s[j] ~ dnorm(%s, 1 / pow(%s, 2))", f, mu, sdf))
  }

  # ---- priors -------------------------------------------------------------------------
  lines <- character(); seen <- character(); ssp_groups <- character(); hier <- character()
  ssp_lines <- function(target, group)
    c(sprintf("  incl_%s ~ dbern(p_incl_%s)", target, group),
      sprintf("  slab_%s ~ dnorm(0, 1 / pow(sigma_slab_%s, 2))", target, group),
      sprintf("  %s <- incl_%s * slab_%s", target, target, target))
  for (i in seq_len(nrow(ld))) {
    r <- ld[i, ]; n <- lnode(r$lhs, r$rhs)
    if (r$fixed != "") {
      add(r$lhs, "=~", r$rhs, fixed = as.numeric(r$fixed), kind = "loading")
      if (!r$zero) lines <- c(lines, sprintf("  %s <- %s", n, r$fixed))
      next
    }
    target <- if (r$label != "") sprintf("lab_%s", safe(r$label)) else n
    pr <- if (r$prior != "") r$prior else if (r$first) dp$loading else dp$cross
    jpr <- pr
    if (pr == "hier") { hier <- union(hier, r$lhs); pr <- sprintf("hier(%s)", r$lhs); jpr <- sprintf("dlnorm(mu_load_%s, 1)", r$lhs) }
    add(r$lhs, "=~", r$rhs, node = target, prior = pr, kind = "loading",
        incl = if (r$prior == "ssp") sprintf("incl_%s", target) else NA_character_)
    if (r$label != "") lines <- c(lines, sprintf("  %s <- %s", n, target))
    if (target %in% seen) next
    seen <- c(seen, target)
    if (r$first && r$prior != "" && r$prior != "ssp" && prior_density(r$prior)$lo < 0)    # measurement loadings are positive (as in RTMB / ELGM)
      jpr <- sprintf("%s T(0,)", r$prior)
    lines <- c(lines, if (r$prior == "") sprintf("  %s ~ %s", target, jpr)
               else if (r$prior == "ssp") { ssp_groups <- union(ssp_groups, "load"); ssp_lines(target, "load") }
               else sprintf("  %s ~ %s", target, jpr))
  }
  betas <- character()
  for (i in seq_len(nrow(rg))) {
    r <- rg[i, ]
    n <- bnode(r$lhs, r$x)
    target <- if (r$label != "") sprintf("labb_%s", safe(r$label)) else n            # shared label: one node
    add(r$lhs, "~", r$x, node = target, prior = if (r$prior == "") dp$beta else r$prior, kind = "beta",
        incl = if (r$prior == "ssp") sprintf("incl_%s", target) else NA_character_)
    if (r$label != "") betas <- c(betas, sprintf("  %s <- %s", n, target))
    if (target %in% seen) next
    seen <- c(seen, target); n <- target
    betas <- c(betas, if (r$prior == "") sprintf("  %s ~ %s", n, dp$beta)
               else if (r$prior == "ssp") { ssp_groups <- union(ssp_groups, "reg"); ssp_lines(n, "reg") }
               else sprintf("  %s ~ %s", n, r$prior))
  }
  items <- character()
  for (v in B) { items <- c(items, sprintf("  d_%s ~ %s", safe(v), dp$intercept_binary))
                 add(v, "~1", "", node = sprintf("d_%s", safe(v)), kind = "intercept", prior = dp$intercept_binary) }
  for (v in C) {
    items <- c(items, sprintf("  xi_%s ~ %s", safe(v), dp$intercept),
               sprintf("  %s_%s ~ %s", if (is_q(fam$ind[[v]])) "sdr" else "sigma", safe(v), dp$resid_sd))
    add(v, "~1", "", node = sprintf("xi_%s", safe(v)), kind = "intercept", prior = dp$intercept)
    add(v, "~~", v, node = sprintf("sigma_%s", safe(v)), kind = sprintf("resid%s", if (is_q(fam$ind[[v]])) paste0("_", fam$ind[[v]]$family) else ""),
        prior = sprintf("sd %s", dp$resid_sd))
  }
  for (v in O) {                                                                     # thresholds: first free, then positive increments
    sv <- safe(v); K <- spec$ord$K[match(v, spec$ord$items)]; r <- ld[ld$rhs == v & !ld$zero, ]; ty <- spec$ord$type[[v]]
    if (ty != "grm") {                                                               # partial credit: free step difficulties
      items <- c(items, sprintf("  for (k in 1:%d) { b_%s[k] ~ %s }", K - 1, sv, dp$threshold))
      for (k in seq_len(K - 1)) add(v, "|", sprintf("t%d", k), node = sprintf("b_%s[%d]", sv, k), kind = "threshold", prior = dp$threshold)
      if (ty == "tppcm") {                                                           # step discriminations a_k = a r_k, r_1 = 1
        items <- c(items, sprintf("  r_%s[1] <- 1", sv))
        if (K > 2) items <- c(items, sprintf("  for (k in 2:%d) { r_%s[k] ~ %s; as_%s[k] <- %s * r_%s[k] }", K - 1, sv, dp$step_ratio, sv, lnode(r$lhs, v), sv))
        for (k in seq_len(K - 2) + 1) add(v, "|a", sprintf("a%d", k), node = sprintf("as_%s[%d]", sv, k), kind = "step_disc", prior = sprintf("loading x %s", dp$step_ratio))
      }
      next
    }
    items <- c(items, sprintf("  b_%s[1] ~ %s", sv, dp$threshold))
    if (K > 2) items <- c(items, sprintf("  for (k in 2:%d) { inc_%s[k] ~ %s; b_%s[k] <- b_%s[k - 1] + inc_%s[k] }", K - 1, sv, dp$threshold_step, sv, sv, sv))
    if (!grepl("[j]", lexpr(r$lhs, v), fixed = TRUE)) items <- c(items, sprintf("  for (k in 1:%d) { cut_%s[k] <- %s * b_%s[k] }", K - 1, sv, lnode(r$lhs, v), sv))
    for (k in seq_len(K - 1)) add(v, "|", sprintf("t%d", k), node = sprintf("b_%s[%d]", sv, k), kind = "threshold",
                                  prior = if (k == 1) dp$threshold else dp$threshold_step)
  }
  modl <- character(); mgroups <- character()
  for (i in seq_len(nrow(mo))) {                                                     # moderation effects
    r <- mo[i, ]; n <- mnode(r); al <- r$kind == "alpha"
    lhs <- sprintf("E(%s)", r$target); rhs <- if (al) paste0(r$mod, ":", r$factor) else r$mod
    kd <- if (al) "alpha" else "delta"; dflt <- if (al) dp$moderation else dp$beta
    g <- sprintf("%s_%s", if (al) "alm" else "dlt", safe(r$mod))
    if (r$type == "none") { add(lhs, "~", rhs, fixed = 0, kind = kd); next }
    add(lhs, "~", rhs, node = n, prior = if (r$type == "ssp") "ssp" else if (r$prior != "") r$prior else dflt, kind = kd,
        incl = if (r$type == "ssp") sprintf("incl_%s", n) else NA_character_, group = if (r$type == "ssp") g else NA_character_)
    if (n %in% seen) next
    seen <- c(seen, n)
    modl <- c(modl, if (r$type == "ssp") { mgroups <- union(mgroups, g); ssp_lines(n, g) }
                    else sprintf("  %s ~ %s", n, if (r$prior != "") r$prior else dflt))
  }
  for (f in fs) for (z in smods(fam$factor[[f]])) {                                 # V(f) ~ z, V(y) ~ z
    pr <- fam$factor[[f]]$scale_prior[z] %||% NA; pr <- if (is.na(pr) || pr == "") dp$moderation else pr
    n <- sprintf("psi_%s_%s", f, safe(z)); modl <- c(modl, sprintf("  %s ~ %s", n, pr)); add(sprintf("V(%s)", f), "~", z, node = n, prior = pr, kind = "psi")
  }
  for (v in C) for (z in smods(fam$ind[[v]])) {
    pr <- fam$ind[[v]]$scale_prior[z] %||% NA; pr <- if (is.na(pr) || pr == "") dp$moderation else pr
    n <- sprintf("kap_%s_%s", safe(v), safe(z)); modl <- c(modl, sprintf("  %s ~ %s", n, pr)); add(sprintf("V(%s)", v), "~", z, node = n, prior = pr, kind = "kappa")
  }
  ssp_groups <- c(ssp_groups, mgroups)
  hyper <- character()
  for (f in fs) {
    ff <- fam$factor[[f]]
    kd <- if (is_q(ff)) paste0("fscale_", ff$family) else "fsd"
    if (is.na(spec$fsd[[f]])) hyper <- c(hyper, sprintf("  sd_%s ~ %s", f, dp$factor_sd))
    fpr <- if (is.na(spec$fsd[[f]])) sprintf("sd %s", dp$factor_sd) else ""
    if (is_q(ff)) add(f, "~~", f, node = sprintf("scale_%s", f), fixed = spec$fsd[[f]], kind = kd, prior = fpr)
    else if (is.na(spec$fsd[[f]])) add(f, "~~", f, node = sprintf("sd_%s", f), kind = kd, prior = fpr)
    else add(f, "~~", f, fixed = spec$fsd[[f]], kind = kd)
  }
  if (isTRUE(spec$cov_free)) {
    n <- sprintf("r_%s_%s", fs[1], fs[2]); hyper <- c(hyper, sprintf("  %s ~ %s", n, dp$cor))
    add(fs[1], "~~", fs[2], node = n, kind = "corr", prior = dp$cor)
  }
  for (f in hier) hyper <- c(hyper, hier_jags(f))
  for (g in ssp_groups) hyper <- c(hyper,
    if (is.null(ssp_p)) sprintf("  p_incl_%s ~ dbeta(1, 1)", g) else sprintf("  p_incl_%s <- %s", g, num(ssp_p)),
    sprintf("  sigma_slab_%s ~ %s", g, dp$slab_sd))

  code <- paste0(
    "# generated by birt ", pkg_version(), "\nmodel {\n  for (j in 1:N) {\n",
    paste(lik, collapse = "\n"), "\n", paste(persons, collapse = "\n"), "\n  }\n\n",
    "  # loadings\n", paste(lines, collapse = "\n"), "\n\n",
    "  # intercepts, thresholds and residual scales\n", paste(items, collapse = "\n"), "\n",
    if (length(betas)) paste0("\n  # latent regression\n", paste(betas, collapse = "\n"), "\n"),
    if (length(modl)) paste0("\n  # moderation (E: shifts and log loadings; V: log SDs)\n", paste(modl, collapse = "\n"), "\n"),
    if (length(hyper)) paste0("\n  # latent variables and hyperparameters\n", paste(hyper, collapse = "\n"), "\n"),
    if (length(qglob)) paste0("\n  # ALD scales from the implied SDs\n", paste(unique(qglob), collapse = "\n"), "\n"), "}\n")

  partab <- do.call(rbind, rows)
  monitor <- unique(c(stats::na.omit(partab$node), stats::na.omit(partab$incl), if (length(O)) sprintf("r_%s", safe(O[spec$ord$type[O] == "tppcm"])),
                      if (is.null(ssp_p)) sprintf("p_incl_%s", ssp_groups), sprintf("sigma_slab_%s", ssp_groups),
                      sprintf("mu_load_%s", hier), fs, sprintf("b_%s", safe(O)),
                      sprintf("as_%s", safe(O[spec$ord$type[O] == "tppcm" & spec$ord$K[match(O, spec$ord$items)] > 2]))))
  list(code = code, partab = partab, monitor = monitor)
}
