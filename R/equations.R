# equations(): the model of a fit written out as equations (Unicode or LaTeX) ---------------
#
# Each model is a list of sections (measurement, structural, priors); a line is
# (lhs, relation, rhs, note). Symbols are rendered by a format object F, so one builder
# gives both the Unicode text (default; for the console) and LaTeX (an aligned block for a
# paper). Symbols: persons j, items i, categories k; lambda measurement loadings, rho
# cross-loadings, d / xi intercepts, sigma residual SDs, beta covariate effects, gamma the path
# between latent variables, zeta latent residuals.

#' Model equations
#'
#' Writes the fitted model as equations, with the parameters as they appear in [summary()],
#' the sets of items they apply to and, for Bayesian fits, the priors. The default is Unicode
#' text for the console; `format = "latex"` gives an `aligned` block for a paper or an
#' appendix (e.g. `writeLines(as.character(equations(fit, "latex")), "model.tex")`).
#' @param object A fitted model (`birt()` with either engine, or `hgrm()`).
#' @param format `"unicode"` or `"latex"`.
#' @param priors Include the priors (Bayesian fits).
#' @param x A `birt_equations` object.
#' @param ... Unused.
#' @return A `birt_equations` object (printed); `as.character()` gives the lines.
#' @examples
#' cat(rtirt_syntax(paste0("y", 1:3), paste0("t", 1:3)))
#' \donttest{
#' if (requireNamespace("RTMB", quietly = TRUE)) {
#' d <- sim_rtirt(N = 200, K = 3)
#' fit <- birt(rtirt_syntax(paste0("y", 1:3), paste0("t", 1:3)), d,
#'             control = rtmb_control(method = "aghq"), progress = FALSE)
#' equations(fit)
#' equations(fit, "latex")
#' }
#' }
#' @export
equations <- function(object, ...) UseMethod("equations")

# ---- rendering ----------------------------------------------------------------------------

eq_fmt <- function(format) {
  u <- format == "unicode"
  gr <- c(alpha = "\u03b1", beta = "\u03b2", gamma = "\u03b3", delta = "\u03b4", epsilon = "\u03b5", zeta = "\u03b6",
          eta = "\u03b7", theta = "\u03b8", kappa = "\u03ba", lambda = "\u03bb", nu = "\u03bd", xi = "\u03be",
          mu = "\u03bc", pi = "\u03c0", rho = "\u03c1", sigma = "\u03c3", tau = "\u03c4", psi = "\u03c8")
  usub <- c(i = "\u1d62", j = "\u2c7c", k = "\u2096", p = "\u209a", `0` = "\u2080", `1` = "\u2081", `2` = "\u2082",
            `3` = "\u2083", `4` = "\u2084", `5` = "\u2085", `6` = "\u2086", `7` = "\u2087", `8` = "\u2088", `9` = "\u2089",
            `-` = "\u208b", `+` = "\u208a", `,` = ",")
  sub1 <- function(s) {
    ch <- strsplit(as.character(s), "")[[1]]
    if (all(ch %in% names(usub))) paste(usub[ch], collapse = "") else paste0("[", s, "]")
  }
  list(u = u,
    g = function(n) if (u) gr[[n]] else paste0("\\", n),
    v = function(n) if (u || nchar(n) == 1) n else sprintf("\\mathrm{%s}", gsub("_", "\\\\_", n)),
    s = function(x, s) if (u) paste0(x, sub1(s)) else sprintf("%s_{%s}", x, s),
    sim = if (u) "\u223c" else "\\sim", eq = "=", dot = if (u) "\u00b7" else "\\,", minus = if (u) "\u2212" else "-",
    sq = function(x) if (x == "1") "1" else if (u) paste0(x, "\u00b2") else paste0(x, "^2"),
    inset = if (u) " \u2208 " else " \\in ",
    N = function(m, v, plus = FALSE) sprintf(if (u) "N%s(%s, %s)" else "\\mathcal{N}%s(%s, %s)", if (!plus) "" else if (u) "\u207a" else "^{+}", m, v),
    fun = function(f, x) if (u) sprintf("%s(%s)", f, x) else sprintf("\\%s(%s)", f, x),
    op = function(f) if (u) f else if (f %in% c("log", "exp")) paste0("\\", f) else sprintf("\\operatorname{%s}", f),
    text = function(t) if (u) t else sprintf("\\text{%s}", gsub("_", "\\\\_", t)),
    set = function(items) {
      it <- if (length(items) > 4) c(items[1:2], if (u) "\u2026" else "\\ldots", items[length(items)]) else items
      if (u) sprintf("{%s}", paste(it, collapse = ", ")) else sprintf("\\{%s\\}", paste(vapply(it, function(x) if (x == "\\ldots") x else sprintf("\\mathrm{%s}", x), ""), collapse = ", "))
    })
}

eq_line <- function(lhs, rel, rhs, note = "") list(lhs = lhs, rel = rel, rhs = rhs, note = note)

eq_object <- function(sections, format) structure(list(sections = Filter(function(s) length(s$lines) > 0, sections), format = format),
                                                  class = "birt_equations")

#' @rdname equations
#' @export
as.character.birt_equations <- function(x, ...) {
  if (x$format == "unicode") {
    out <- character()
    for (s in x$sections) {
      out <- c(out, if (length(out)) "", paste0(s$title, ":"))
      L <- s$lines; wl <- max(nchar(vapply(L, `[[`, "", "lhs"), type = "width"))
      for (l in L) {
        pad <- strrep(" ", wl - nchar(l$lhs, type = "width"))
        out <- c(out, sub("\\s+$", "", sprintf("  %s%s %s %s%s", pad, l$lhs, l$rel, l$rhs, if (nzchar(l$note)) paste0("    ", l$note) else "")))
      }
    }
    return(out)
  }
  body <- unlist(lapply(x$sections, function(s)
    c(sprintf("  & \\text{%s:} \\\\", s$title),
      vapply(s$lines, function(l) sprintf("  %s &%s %s%s \\\\", l$lhs, if (nzchar(l$rel)) l$rel else "\\quad", l$rhs,
                                           if (nzchar(l$note)) sprintf(" \\quad %s", l$note) else ""), ""))))
  body[length(body)] <- sub(" \\\\\\\\$", "", body[length(body)])
  c("\\begin{aligned}", body, "\\end{aligned}")
}

#' @rdname equations
#' @export
print.birt_equations <- function(x, ...) { cat(as.character(x), sep = "\n"); invisible(x) }

# JAGS prior string -> math notation (precision -> SD); falls back to the string itself
prior_math <- function(txt, F) {
  txt <- trimws(txt)
  if (!nzchar(txt)) return("")
  if (txt == "ssp") return(F$text("spike-and-slab"))
  if (grepl("^hier\\(", txt))                                   # hierarchical lognormal loading prior
    return(sprintf(if (F$u) "logN(%s, 1)" else "\\log\\mathcal{N}(%s, 1)", F$s(F$g("mu"), F$v(sub("^hier\\((.*)\\)$", "\\1", txt)))))
  sd_pr <- grepl("^sd ", txt); txt <- sub("^sd ", "", txt)
  m <- regmatches(txt, regexec("^(d[a-z]+)\\(([^()]*)\\)\\s*(T\\(([^()]*)\\))?$", txt))[[1]]
  if (!length(m)) return(if (F$u) txt else sprintf("\\texttt{%s}", txt))
  a <- trimws(strsplit(m[3], ",")[[1]]); num <- function(z) tryCatch(eval(str2lang(z), baseenv()), error = function(e) NA)
  fmt <- function(v) format(signif(v, 3))
  tr <- m[5]; pos <- grepl("^\\s*0\\s*,\\s*$", tr)
  sdv <- function(tau) { s <- 1 / sqrt(num(tau)); if (is.na(s)) tau else fmt(s) }
  out <- switch(m[2],
    dnorm = F$N(a[1], F$sq(sdv(a[2])), plus = pos),
    dlnorm = sprintf(if (F$u) "logN(%s, %s)" else "\\log\\mathcal{N}(%s, %s)", sub("^-", F$minus, a[1]), F$sq(sdv(a[2]))),
    dt = sprintf(if (F$u) "t%s%s(%s, %s)" else "t^{%s}_{%s}(%s, %s)", if (pos) (if (F$u) "\u207a" else "+") else "",
                 if (F$u) F$s("", a[3]) else a[3], a[1], sdv(a[2])),
    dunif = sprintf(if (F$u) "U(%s, %s)" else "\\mathcal{U}(%s, %s)", a[1], a[2]),
    dbeta = sprintf(if (F$u) "Beta(%s, %s)" else "\\mathrm{Beta}(%s, %s)", a[1], a[2]),
    dgamma = sprintf(if (F$u) "Gamma(%s, %s)" else "\\mathrm{Gamma}(%s, %s)", a[1], a[2]),
    dexp = sprintf(if (F$u) "Exp(%s)" else "\\mathrm{Exp}(%s)", a[1]),
    if (F$u) txt else sprintf("\\texttt{%s}", txt))
  if (nzchar(tr) && !pos && m[2] != "dt") out <- paste0(out, sprintf(if (F$u) " T(%s)" else "\\,T(%s)", tr))
  if (m[2] == "dt" && nzchar(tr) && !pos) out <- paste0(out, sprintf(" T(%s)", tr))
  attr(out, "sd") <- sd_pr
  out
}

# items -> a common letter (y1..y4 -> y), or a fallback
item_letter <- function(items, fallback) {
  p <- unique(sub("[0-9_.]+$", "", items))
  if (length(p) == 1 && grepl("^[A-Za-z]{1,3}$", p)) p else fallback
}

# ---- lavaan-syntax models -------------------------------------------------------------------

eq_lavaan <- function(sp, e, bayes, format, priors, p_ald = NULL, marginal = FALSE) {
  F <- eq_fmt(format); fs <- sp$factors
  I <- F$s(F$v("x"), "j")
  fsym <- function(f) F$s(F$v(f), "j")
  ld <- sp$ld[!sp$ld$zero, ]
  # measurement: group indicators with the same structure
  mo <- sp$mod %||% mod_empty()
  has <- function(v, kind) any(mo$target == v & mo$kind == kind & mo$type != "none")
  key <- vapply(sp$ovs, function(v) {
    r <- ld[ld$rhs == v, ]
    fam <- if (v %in% sp$binary) "bin" else if (v %in% (sp$ordinal %||% character())) paste0("ord", sp$ord$type[[v]]) else { f <- sp$family$ind[[v]]; paste(f$family, f$p, paste(f$scale, collapse = "+"), paste(f$skew, collapse = "+")) }
    paste(fam, has(v, "alpha"), has(v, "beta"), paste(sprintf("%s:%s%s", r$lhs, ifelse(r$fixed != "", r$fixed, "free"), ifelse(r$first, "L", "C")), collapse = ","))
  }, "")
  meas <- list(); modl <- list()
  for (k in unique(key)) {
    vs <- sp$ovs[key == k]; bin <- vs[1] %in% sp$binary; ord <- vs[1] %in% (sp$ordinal %||% character()); r <- ld[ld$rhs == vs[1], ]
    al <- has(vs[1], "alpha"); be <- has(vs[1], "beta")
    yl <- item_letter(vs, if (bin || ord) "y" else "w")
    note <- sprintf("i%s%s", F$inset, F$set(vs))
    dl <- F$s(F$g("delta"), "ij")
    if (ord) {                                                    # grm: P(y > k); gpcm / tppcm: adjacent categories
      ty <- sp$ord$type[[vs[1]]]
      lam <- F$s(F$g("lambda"), paste0(if (al) "ij" else "i", if (ty == "tppcm") "k" else ""))
      lhs <- if (ty == "grm") sprintf("%s P(%s > k)", F$op("logit"), F$s(F$v(yl), "ij"))
             else sprintf("%s P(%s = k + 1) / P(%s = k)", F$op("log"), F$s(F$v(yl), "ij"), F$s(F$v(yl), "ij"))
      meas[[length(meas) + 1]] <- eq_line(lhs, F$eq,
        sprintf("%s(%s %s %s%s)", lam, fsym(r$lhs[1]), F$minus, F$s("b", "ik"), if (be) paste0(" ", F$minus, " ", dl) else ""),
        sprintf("%s, k = 1, %s, K%s %s 1%s", note, if (F$u) "\u2026" else "\\ldots", if (F$u) "\u1d62" else "_i", F$minus,
                switch(ty, grm = "", gpcm = paste0("; ", F$text("generalized partial credit")),
                       tppcm = paste0("; ", F$text("two-parameter partial credit"), ", ", F$s(F$g("lambda"), "i1"), " = ", F$s(F$g("lambda"), "i")))))
      next
    }
    terms <- character()
    for (q in seq_len(nrow(r))) {
      f <- r$lhs[q]; fx <- r$fixed[q]
      sym <- if (r$first[q]) F$s(F$g("lambda"), if (al) "ij" else "i") else F$s(F$g("rho"), "i")
      t <- if (fx != "") {
        v <- as.numeric(fx); if (v == 1) fsym(f) else if (v == -1) paste0(F$minus, " ", fsym(f)) else paste0(format(v), F$dot, fsym(f))
      } else paste0(sym, F$dot, fsym(f))
      terms <- c(terms, t)
    }
    if (be && bin) terms[1] <- sprintf("%s(%s %s %s)", F$s(F$g("lambda"), if (al) "ij" else "i"), fsym(r$lhs[1]), F$minus, dl)
    lin <- paste(c(F$s(if (bin) "d" else F$g("xi"), "i"), terms, if (be && !bin) dl), collapse = " + ")
    lin <- gsub(paste0("\\+ ", F$minus), F$minus, lin, fixed = FALSE)
    lin <- gsub(paste0("+ ", F$minus), F$minus, lin, fixed = TRUE)
    if (bin) meas[[length(meas) + 1]] <- eq_line(sprintf("%s P(%s = 1)", F$op("logit"), F$s(F$v(yl), "ij")), F$eq, lin, note)
    else {
      fam <- sp$family$ind[[vs[1]]]; eps <- F$s(F$g("epsilon"), "ij")
      meas[[length(meas) + 1]] <- eq_line(F$s(F$v(yl), "ij"), F$eq, paste(lin, "+", eps), note)
      sg <- if (length(fam$scale)) F$s(F$g("sigma"), "ij") else F$s(F$g("sigma"), "i")
      dist <- switch(fam$family,
        normal = F$N("0", F$sq(sg)),
        ald = sprintf(if (F$u) "ALD%s(0, %s)" else "\\mathrm{ALD}_{%s}(0, %s)", if (F$u) paste0("[p = ", fam$p, "]") else paste0("p=", fam$p), sg),
        shash = sprintf(if (F$u) "SHASH(0, %s, %s, %s)" else "\\mathrm{SHASH}(0, %s, %s, %s)", sg,
                        if (length(fam$skew)) F$s(F$g("kappa"), "ij") else F$s(F$g("kappa"), "i"), F$s(F$g("delta"), "i")))
      meas[[length(meas) + 1]] <- eq_line(eps, F$sim, dist, if (fam$family == "shash") F$text("skew, tail weight (below 1: heavier tails than normal)") else "")
      if (length(fam$scale)) meas[[length(meas) + 1]] <- eq_line(paste(F$op("log"), sg), F$eq,
        paste(c(paste(F$op("log"), F$s(F$g("sigma"), "i")), vapply(fam$scale, function(z) paste0(F$s(F$g("psi"), paste0("i")), F$dot, F$s(F$v(z), "j")), "")), collapse = " + "), "")
      if (length(fam$skew)) meas[[length(meas) + 1]] <- eq_line(F$s(F$g("kappa"), "ij"), F$eq,
        paste(c(F$s(F$g("kappa"), "i"), vapply(fam$skew, function(z) paste0(F$s(F$g("eta"), "i"), F$dot, F$s(F$v(z), "j")), "")), collapse = " + "), "")
    }
  }
  # structural
  st <- list()
  for (f in fs) {
    r <- sp$reg[sp$reg$lhs == f, ]
    terms <- vapply(seq_len(nrow(r)), function(q) {
      x <- r$x[q]
      if (r$is_factor[q]) paste0(F$g("gamma"), F$dot, fsym(x)) else paste0(F$s(F$g("beta"), if (F$u) x else sprintf("\\mathrm{%s}", x)), F$dot, F$s(F$v(x), "j"))
    }, "")
    fam <- sp$family$factor[[f]]; z <- F$s(F$g("zeta"), "j")
    sdf <- sp$fsd[[f]]
    sv <- if (is.na(sdf)) F$sq(F$s(F$g("sigma"), F$v(f))) else format(sdf^2)
    sg <- if (length(fam$scale)) F$s(F$g("sigma"), "j") else if (is.na(sdf)) F$s(F$g("sigma"), F$v(f)) else format(sdf)
    dist <- switch(fam$family,
      normal = F$N("0", if (length(fam$scale)) F$sq(sg) else sv),
      ald = sprintf(if (F$u) "ALD[p = %s](0, %s)" else "\\mathrm{ALD}_{p=%s}(0, %s)", fam$p, if (is.na(sdf)) sg else F$text("scale with variance 1")),
      shash = sprintf(if (F$u) "SHASH(0, %s, %s, %s)" else "\\mathrm{SHASH}(0, %s, %s, %s)", sg, F$g("kappa"), F$g("delta")))
    if (length(terms)) {
      st[[length(st) + 1]] <- eq_line(fsym(f), F$eq, paste(c(terms, z), collapse = " + "), "")
      st[[length(st) + 1]] <- eq_line(z, F$sim, dist, "")
    } else st[[length(st) + 1]] <- eq_line(fsym(f), F$sim, sub("^(\\\\mathcal\\{N\\}|N)\\(0", "\\1(0", dist), "")
    if (length(fam$scale)) st[[length(st) + 1]] <- eq_line(paste(F$op("log"), sg), F$eq,
      paste(c(paste(F$op("log"), if (is.na(sdf)) F$s(F$g("sigma"), F$v(f)) else format(sdf)), vapply(fam$scale, function(z) paste0(F$s(F$g("psi"), F$v(z)), F$dot, F$s(F$v(z), "j")), "")), collapse = " + "), "")
  }
  # moderation of loadings and shifts: the equations, then which indicators carry which effect
  act <- mo[mo$type != "none", ]
  zs <- function(rr) paste(vapply(seq_len(nrow(rr)), function(q) paste0(F$s(F$g(if (rr$kind[q] == "alpha") "alpha" else "beta"), if (length(unique(act$mod)) > 1) paste0("i,", rr$mod[q]) else "i"), F$dot, F$s(F$v(rr$mod[q]), "j")), ""), collapse = " + ")
  if (any(act$kind == "alpha")) modl[[length(modl) + 1]] <- eq_line(paste(F$op("log"), F$s(F$g("lambda"), "ij")), F$eq,
    paste(paste(F$op("log"), F$s(F$g("lambda"), "i")), zs(act[act$kind == "alpha" & !duplicated(paste(act$kind, act$mod)), ]), sep = " + "), F$text("E(y) ~ z:f"))
  if (any(act$kind == "beta")) modl[[length(modl) + 1]] <- eq_line(F$s(F$g("delta"), "ij"), F$eq,
    zs(act[act$kind == "beta" & !duplicated(paste(act$kind, act$mod)), ]), F$text("E(y) ~ z: on the latent scale for categorical indicators"))
  for (kd in c("beta", "alpha")) for (z in unique(mo$mod[mo$kind == kd])) for (ty in c("ssp", "free", "common", "none")) {
    j <- mo$target[mo$kind == kd & mo$mod == z & mo$type == ty]; if (!length(j)) next
    what <- switch(ty, ssp = F$text("spike-and-slab"), free = F$text("free"), none = paste0("0 ", F$text("(anchor)")),
                   common = F$text(sprintf("equal, label %s", paste(unique(mo$label[mo$kind == kd & mo$mod == z & mo$type == ty]), collapse = ", "))))
    modl[[length(modl) + 1]] <- eq_line(F$s(F$g(if (kd == "alpha") "alpha" else "beta"), paste0("i,", z)), ":", what, sprintf("i%s%s", F$inset, F$set(j)))
  }
  if (isTRUE(sp$cov_free)) st[[length(st) + 1]] <- eq_line(sprintf("%s(%s, %s)", F$op("Corr"), fsym(fs[1]), fsym(fs[2])), F$eq, F$v("r"), "")
  secs <- list(list(title = "Measurement", lines = meas), list(title = "Latent variables", lines = st), list(title = "Moderation", lines = modl),
               list(title = "Censoring (censor =)", lines = eq_censor(sp$censor, F)))
  if (priors && bayes) secs[[length(secs) + 1]] <- list(title = "Priors", lines = eq_priors_lavaan(sp, e, F))
  secs[[length(secs) + 1]] <- eq_likelihood(F, bayes, marginal)
  secs
}

# censored continuous indicators: the likelihood contribution of a value beyond its limit
eq_censor <- function(cs, F) {
  if (is.null(cs)) return(list())
  y <- F$s(F$v(item_letter(cs$item, "w")), "ij"); fy <- function(c) sprintf(if (F$u) "F\u1d62(%s | %s)" else "F_i(%s \\mid %s)", c, F$s(F$v("f"), "j"))
  lim <- function(v, raw, lg) if (lg) paste(format(signif(v, 4)), F$text(sprintf("(log %s)", format(raw)))) else format(signif(v, 4))
  out <- list()
  for (side in c("upper", "lower")) {
    r <- cs[is.finite(cs[[side]]), ]; if (!nrow(r)) next
    key <- sprintf("%s|%s", r[[side]], r$log)
    for (k in unique(key)) {
      i <- which(key == k); c0 <- lim(r[[side]][i[1]], r[[paste0(side, "_raw")]][i[1]], r$log[i[1]])
      cc <- F$s(F$v("c"), "i")
      out[[length(out) + 1]] <- eq_line(sprintf("%s %s %s", y, if (side == "upper") (if (F$u) "\u2265" else "\\ge") else (if (F$u) "\u2264" else "\\le"), cc), ":",
        if (side == "upper") sprintf("P(%s %s %s) = 1 %s %s", y, if (F$u) "\u2265" else "\\ge", cc, F$minus, fy(cc))
        else sprintf("P(%s %s %s) = %s", y, if (F$u) "\u2264" else "\\le", cc, fy(cc)),
        sprintf("i%s%s; %s = %s; %s", F$inset, F$set(r$item[i]), cc, c0, F$text(if (side == "upper") "right censored: only 'at least c' is known" else "left censored: only 'at most c' is known")))
    }
  }
  out
}

# which likelihood a fit uses: conditional (persons sampled, JAGS) or marginal (persons integrated out, RTMB)
eq_likelihood <- function(F, bayes, marginal) {
  th <- F$s(F$g("theta"), "j"); ps <- F$g("psi")
  cond <- if (F$u) sprintf("p(%s, %s | y) \u221d \u220f\u2c7c p(y\u2c7c | %s, %s) p(%s | %s) p(%s)", ps, "\u03b8", th, ps, th, ps, ps)
          else sprintf("p(\\psi, \\theta \\mid y) \\propto \\prod_j p(y_j \\mid %s, \\psi)\\, p(%s \\mid \\psi)\\, p(\\psi)", th, th)
  marg <- if (F$u) sprintf("L(%s) = \u220f\u2c7c \u222b p(y\u2c7c | %s, %s) p(%s | %s) d%s", ps, th, ps, th, ps, th)
          else sprintf("L(\\psi) = \\prod_j \\int p(y_j \\mid %s, \\psi)\\, p(%s \\mid \\psi)\\, d%s", th, th, th)
  lines <- if (!marginal) list(eq_line(F$text("conditional likelihood"), "", cond, F$text("persons sampled with the parameters (JAGS); psi: all other parameters")))
           else list(eq_line(F$text("marginal likelihood"), "", marg, F$text(if (bayes) "persons integrated out (AGHQ); posterior of psi by ELGM: p(psi | y) proportional to L(psi) p(psi)"
                                                                            else "persons integrated out by adaptive Gauss-Hermite quadrature; maximum likelihood")))
  list(title = "Likelihood", lines = lines)
}

eq_priors_lavaan <- function(sp, e, F) {
  sym_of <- function(r) {
    if (r$op == "=~") { first <- sp$ld$first[sp$ld$lhs == r$lhs & sp$ld$rhs == r$rhs][1]; F$s(F$g(if (isTRUE(first)) "lambda" else "rho"), "i") }
    else if (r$op == "~1") F$s(if (r$lhs %in% sp$binary) "d" else F$g("xi"), "i")
    else if (r$op == "~~" && r$lhs == r$rhs && r$lhs %in% sp$ovs) F$s(F$g("sigma"), "i")
    else if (r$op == "~~" && r$lhs == r$rhs) F$s(F$g("sigma"), F$v(r$lhs))
    else if (r$op == "~~") F$v("r")
    else if (r$op == "~" && grepl("^E\\(", r$lhs)) F$s(F$g(if (grepl(":", r$rhs)) "alpha" else "beta"), paste0("i,", sub(":.*", "", r$rhs)))
    else if (r$op == "~" && grepl("^V\\(", r$lhs)) F$s(F$g("psi"), if (sub("^V\\((.*)\\)$", "\\1", r$lhs) %in% sp$factors) F$v(r$rhs) else paste0("i,", r$rhs))
    else if (r$op == "~") if (r$rhs %in% sp$factors) F$g("gamma") else F$s(F$g("beta"), if (F$u) r$rhs else sprintf("\\mathrm{%s}", r$rhs))
    else if (r$op == "|" && r$prior == (sp$dp %||% dpriors())$threshold) if (r$rhs == "t1") F$s("b", "i1") else F$s("b", "ik")
    else if (r$op == "|") paste0(F$s("b", "ik"), " ", F$minus, " ", F$s("b", "i,k-1"))
    else if (r$op == "|a") paste0(F$s(F$g("lambda"), "ik"), " / ", F$s(F$g("lambda"), "i1"))
    else NA_character_
  }
  e <- e[!is.na(e$prior) & nzchar(e$prior) & !(e$fixed %||% FALSE), ]
  if (!nrow(e)) return(list())
  sy <- vapply(seq_len(nrow(e)), function(i) sym_of(e[i, ]), "")
  it <- ifelse(e$op == "=~", e$rhs, ifelse((e$op %in% c("~1", "~~") & e$lhs == e$rhs & e$lhs %in% sp$ovs) | e$op == "|", e$lhs,
               ifelse(e$op == "~" & grepl("^[EV]\\(", e$lhs), sub("^[EV]\\((.*)\\)$", "\\1", e$lhs), NA)))
  key <- paste(sy, e$prior)
  out <- list()
  for (k in unique(key)) {
    i <- which(key == k); pr <- e$prior[i[1]]; pm <- prior_math(pr, F)
    note <- if (any(!is.na(it[i]))) sprintf("i%s%s", F$inset, F$set(it[i])) else ""
    if (pr == "ssp") {
      pm <- sprintf(if (F$u) "%s%s%s, %s %s Bern(%s), %s %s N(0, %s)" else "%s\\,%s%s, \\ %s %s \\mathrm{Bern}(%s), \\ %s %s \\mathcal{N}(0, %s)",
                    F$s(F$g("delta"), "i"), "", F$s("s", "i"), F$s(F$g("delta"), "i"), F$sim, F$g("pi"), F$s("s", "i"), F$sim,
                    F$sq(F$s(F$g("sigma"), F$v("slab"))))
      out[[length(out) + 1]] <- eq_line(sy[i[1]], F$eq, pm, note)
    } else out[[length(out) + 1]] <- eq_line(sy[i[1]], F$sim, pm, note)
  }
  for (f in unique(sub("^hier\\((.*)\\)$", "\\1", grep("^hier\\(", e$prior, value = TRUE)))) {
    out[[length(out) + 1]] <- eq_line(F$s(F$g("mu"), F$v(f)), F$sim, prior_math(sprintf("dnorm(0, %s)", num(1 / HIER_V)), F), F$text("hierarchical loading prior"))
  }
  if (any(e$prior == "ssp")) {
    dp <- sp$dp %||% dpriors()
    out[[length(out) + 1]] <- eq_line(F$g("pi"), if (is.null(sp$ssp_p)) F$sim else F$eq,
                                      if (is.null(sp$ssp_p)) prior_math("dbeta(1, 1)", F) else format(sp$ssp_p), F$text("shared within loadings / within regressions"))
    out[[length(out) + 1]] <- eq_line(F$s(F$g("sigma"), F$v("slab")), F$sim, prior_math(dp$slab_sd, F), "")
  }
  out
}

#' @rdname equations
#' @export
equations.birt <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format)
  eq_object(eq_lavaan(object$spec, object$estimates, TRUE, format, priors), format)
}

#' @rdname equations
#' @export
equations.birt_rtmb <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format)
  bayes <- !is.null(object$rtmb$elgm)
  eq_object(eq_lavaan(object$spec, object$estimates, bayes, format, priors, marginal = TRUE), format)
}

#' @rdname equations
#' @export
equations.birt_quantile <- function(object, format = c("unicode", "latex"), priors = TRUE, ...) {
  format <- match.arg(format)
  eq <- equations(object$fits[[1]], format, priors)
  F <- eq_fmt(format)
  eq$sections <- c(eq$sections, list(list(title = "Quantiles", lines = list(eq_line(F$v("p"), trimws(F$inset), F$set(format(object$p)),
                                                                                    F$text(sprintf("one fit per quantile of %s", object$target)))))))
  eq
}
