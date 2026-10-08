# algorithm(): the computation behind a fit, step by step (Unicode or LaTeX) ----------------
#
# Same renderer as equations(): sections of lines (lhs, relation, rhs, note). Each section is a
# numbered step of what the code does for this fit (JAGS: the nodes and the samplers JAGS
# assigns; RTMB: the Laplace / AGHQ marginal likelihood; ELGM: the inner and outer quadrature),
# with the R / RTMB function that does it in the note. Non-ASCII symbols are written as
# \u{...} escapes; F$t2(unicode, latex) picks the version of the format.

#' The estimation algorithm of a fit
#'
#' Writes out how a fit was computed, as numbered steps with their formulas: for JAGS fits
#' the full conditionals and the sampler JAGS assigned to every node (from
#' `rjags::list.samplers()` on the compiled model) and the MCMC settings; for RTMB fits the
#' marginal likelihood by the Laplace approximation or adaptive Gauss-Hermite quadrature
#' (centres, nodes, optimization rounds, standard errors); for `method = "elgm"` the inner
#' (per-person) and outer (parameter) quadrature, the node weights, the log marginal likelihood
#' and how the posterior summaries are formed. The numbers (nodes, rounds, hyperparameters)
#' are those of the fit. The model itself is given by [equations()].
#'
#' The same generic is exported by birtRcpp (Gibbs full conditionals and EM steps of its
#' models); either package's `algorithm()` dispatches to the other's fits.
#' @param object A fit of [birt()] (either engine) or [hgrm()].
#' @param format `"unicode"` or `"latex"`.
#' @param samplers JAGS fits: compile the model (without sampling) to list the samplers JAGS
#'   assigns (needs rjags; a second or so).
#' @param ... Unused.
#' @return A `birt_algorithm` object (printed like [equations()]; `as.character()` gives the
#'   lines).
#' @examples
#' \donttest{
#' if (requireNamespace("RTMB", quietly = TRUE)) {
#'   d <- sim_rtirt(N = 200, K = 3)
#'   m <- rtirt_syntax(paste0("y", 1:3), paste0("t", 1:3))
#'   fit <- birt(m, d, control = rtmb_control(method = "aghq"), progress = FALSE)
#'   algorithm(fit)
#'   algorithm(fit, "latex")
#' }
#' }
#' @export
algorithm <- function(object, ...) UseMethod("algorithm")

#' @rdname algorithm
#' @export
algorithm.default <- function(object, ...) {
  if (inherits(object, c("rtirt", "rtirt_qset")) && isNamespaceLoaded("birtRcpp") &&         # a birtRcpp model: its own function
      "algorithm" %in% getNamespaceExports("birtRcpp"))
    return(getExportedValue("birtRcpp", "algorithm")(object, ...))
  stop("algorithm() needs a fit of birt() or hgrm()", if (inherits(object, c("rtirt", "rtirt_qset"))) " (or a birtRcpp version with algorithm())", call. = FALSE)
}

# ---- helpers ------------------------------------------------------------------------------

alg_fmt <- function(format) {
  F <- eq_fmt(format)
  esc <- function(t) { t <- gsub("([&%#$_{}])", "\\\\\\1", t); t <- gsub("~", "\\\\textasciitilde{}", t); gsub("\\^", "\\\\textasciicircum{}", t) }
  F$t2 <- function(u, l) if (F$u) u else l
  F$tt <- function(t) if (F$u) t else sprintf("\\text{%s}", esc(t))
  F$code <- function(t) if (F$u) t else sprintf("\\texttt{%s}", esc(t))
  F
}

alg_object <- function(sections, format) {
  sections <- lapply(sections, function(s) { s$lines <- Filter(Negate(is.null), s$lines); s })
  sections <- Filter(function(s) length(s$lines) > 0, sections)
  for (i in seq_along(sections)) sections[[i]]$title <- sprintf("%d. %s", i, sections[[i]]$title)
  structure(list(sections = sections, format = format), class = c("birt_algorithm", "birt_equations"))
}

# a text line (no formula); the note names the function
al_txt <- function(F, txt, note = "") eq_line("", "", F$tt(txt), if (nzchar(note)) F$code(note) else "")
al_eq <- function(F, lhs, rel, rhs, note = "") eq_line(lhs, rel, rhs, if (nzchar(note)) F$code(note) else "")

# ---- JAGS ---------------------------------------------------------------------------------

# samplers JAGS assigns to the compiled model (glm and dic modules loaded, as R2jags does)
jags_samplers <- function(code, data) {
  if (!requireNamespace("rjags", quietly = TRUE)) return("the rjags package is not installed")
  before <- rjags::list.modules()
  for (m in c("glm", "dic")) rjags::load.module(m, quiet = TRUE)
  on.exit(for (m in setdiff(c("glm", "dic"), before)) try(rjags::unload.module(m, quiet = TRUE), silent = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".jags"); writeLines(code, tf); on.exit(unlink(tf), add = TRUE)
  jm <- tryCatch(suppressWarnings(rjags::jags.model(tf, data = data, n.chains = 1, n.adapt = 0, quiet = TRUE)),
                 error = function(e) conditionMessage(e))
  if (is.character(jm)) return(paste("the model could not be compiled:", jm))
  rjags::list.samplers(jm)
}

# node names -> short labels: indexed nodes as base[i:j], names differing only in their last
# _suffix as stem_* (n)
node_labels <- function(nodes) {
  base <- sub("\\[.*$", "", nodes)
  idx <- suppressWarnings(as.integer(sub("^[^[]*\\[([0-9]+)\\]$", "\\1", nodes)))
  out <- character()
  for (b in unique(base[grepl("\\[", nodes)])) {
    i <- idx[base == b]
    out <- c(out, if (anyNA(i)) sprintf("%s[ ] (%d)", b, length(i)) else if (length(i) == 1) sprintf("%s[%d]", b, i) else {
      i <- sort(i); if (all(diff(i) == 1)) sprintf("%s[%d:%d]", b, i[1], i[length(i)]) else sprintf("%s[ ] (%d)", b, length(i))
    })
  }
  sc <- unique(base[!grepl("\\[", nodes)])
  stem <- ifelse(grepl("_", sc), sub("_[^_]+$", "_*", sc), sc)
  for (s in unique(stem)) {
    m <- sc[stem == s]
    out <- c(out, if (length(m) > 1) sprintf("%s (%d)", s, length(m)) else m)
  }
  out
}

alg_samplers <- function(F, s) {
  if (is.character(s)) return(list(al_txt(F, "JAGS assigns a sampler to every node when it compiles the model; not listed here:"), al_txt(F, s)))
  if (!length(s)) return(list(al_txt(F, "no stochastic nodes to update")))
  nm <- names(s); blk <- lengths(s) > 1
  out <- list()
  for (k in which(blk)) out[[length(out) + 1]] <- eq_line(F$code(nm[k]), ":", F$tt(paste(node_labels(s[[k]]), collapse = ", ")),
                                                          F$tt(sprintf("one block of %d nodes", length(s[[k]]))))
  for (n in unique(nm[!blk])) {
    nodes <- unlist(s[!blk & nm == n])
    out[[length(out) + 1]] <- eq_line(F$code(n), ":", F$tt(paste(node_labels(nodes), collapse = ", ")),
                                      F$tt(sprintf("%d node%s, one at a time", length(nodes), if (length(nodes) > 1) "s" else "")))
  }
  out
}

alg_jags <- function(F, code, data, mc, samplers, ald = FALSE, extra = NULL) {
  burn <- max(mc$n_burn, mc$n_iter / 2); kept <- (mc$n_iter - mc$n_burn) %/% mc$thin
  pa <- F$t2("p(v | pa(v)) \u{220f}_{c \u{2208} ch(v)} p(c | pa(c))", "p(v \\mid \\mathrm{pa}(v)) \\prod_{c \\in \\mathrm{ch}(v)} p(c \\mid \\mathrm{pa}(c))")
  secs <- list(
    list(title = "Joint posterior", lines = list(
      al_eq(F, F$t2("p(\u{03d1}, f | y)", "p(\\vartheta, f \\mid y)"), F$t2("\u{221d}", "\\propto"),
            F$t2("\u{220f}\u{2c7c} p(y\u{2c7c} | f\u{2c7c}, \u{03d1}) p(f\u{2c7c} | \u{03d1}) \u{00b7} p(\u{03d1})",
                 "\\prod_j p(y_j \\mid f_j, \\vartheta)\\, p(f_j \\mid \\vartheta)\\, p(\\vartheta)"), "jags_code(fit); equations(fit)"),
      al_txt(F, paste0("f: the person latent variables (sampled, not integrated out); ", F$t2("\u{03d1}", "theta"), ": all other parameters")),
      if (ald) al_eq(F, F$t2("x\u{2c7c} | e\u{2c7c}", "x_j \\mid e_j"), F$sim,
                     F$t2("N(\u{03bc}\u{2c7c} + k\u{2081} e\u{2c7c}, k\u{2082} \u{03c3} e\u{2c7c}),  e\u{2c7c} \u{223c} Exp(mean \u{03c3})",
                          "\\mathcal{N}(\\mu_j + k_1 e_j, k_2 \\sigma e_j),\\ e_j \\sim \\mathrm{Exp}(\\text{mean } \\sigma)"),
                     "x: a variable with an ald() family; ALD as a normal-exponential mixture (Kozumi & Kobayashi, 2011); e sampled too"))),
    list(title = "Gibbs sweep (one iteration)", lines = list(
      al_eq(F, F$t2("v | rest", "v \\mid \\text{rest}"), F$t2("\u{221d}", "\\propto"), pa, "full conditional of node (or block) v; ch(v): its children"),
      al_txt(F, "every stochastic node is updated once per iteration from its full conditional, by the sampler JAGS chose for it:"))),
    list(title = "Samplers assigned by JAGS (glm and dic modules loaded, as R2jags does)",
         lines = if (isTRUE(samplers)) alg_samplers(F, jags_samplers(code, data)) else list(al_txt(F, "not listed (samplers = FALSE)"))),
    list(title = "MCMC settings", lines = list(
      al_txt(F, sprintf("%d chains, each in its own JAGS process (R2jags::jags.parallel), JAGS seed %s", mc$n_chains, format(mc$seed)), "R2jags::jags.parallel"),
      al_txt(F, sprintf("per chain: 100 adaptation iterations, %s burn-in iterations (R2jags: max(n_burn, n_iter / 2)), then %d iterations kept with thinning %d (%d draws)",
                        format(burn), mc$n_iter - mc$n_burn, mc$thin, kept)),
      al_txt(F, "summaries: posterior mean, SD and 2.5% / 97.5% quantiles; R-hat (coda::gelman.diag), effective sample size (coda::effectiveSize)"))))
  c(secs, extra)
}

alg_fit_indices <- function(F) list(title = "Fit indices (after sampling)", lines = list(
  al_eq(F, F$t2("log p(y\u{2c7c} | \u{03d1}\u{207d}\u{02e2}\u{207e})", "\\log p(y_j \\mid \\vartheta^{(s)})"), F$eq,
        F$t2("log \u{222b} p(y\u{2c7c} | f, \u{03d1}\u{207d}\u{02e2}\u{207e}) p(f | \u{03d1}\u{207d}\u{02e2}\u{207e}) df",
             "\\log \\int p(y_j \\mid f, \\vartheta^{(s)})\\, p(f \\mid \\vartheta^{(s)})\\, df"),
        "add_fit_indices(): Gauss-Hermite (Q = 41) / Gauss-Laguerre over the latents, up to 1000 draws"),
  al_txt(F, "marginal DIC, WAIC and PSIS-LOO (loo package) from the pointwise marginal log-likelihood")))

#' @rdname algorithm
#' @export
algorithm.birt <- function(object, format = c("unicode", "latex"), samplers = TRUE, ...) {
  format <- match.arg(format); F <- alg_fmt(format); sp <- object$spec
  ald <- any(vapply(c(sp$family$factor, sp$family$ind), function(f) identical(f$family, "ald"), TRUE))
  extra <- if (!is.null(object$fit_indices)) list(alg_fit_indices(F))
  alg_object(alg_jags(F, object$code, object$data, object$mcmc, samplers, ald, extra), format)
}


#' @rdname algorithm
#' @export
algorithm.birt_quantile <- function(object, format = c("unicode", "latex"), samplers = TRUE, ...) {
  format <- match.arg(format); F <- alg_fmt(format)
  a <- algorithm(object$fits[[1]], format, samplers)
  a$sections <- c(a$sections, list(list(title = sprintf("%d. Quantiles", length(a$sections) + 1), lines = list(
    al_eq(F, F$v("p"), trimws(F$inset), F$set(format(object$p)),
          sprintf("one independent JAGS fit per quantile of %s, each as above (shown: p = %s)", object$target, format(object$p[1])))))))
  a
}

# ---- RTMB: maximum likelihood (Laplace / AGHQ) and ELGM -----------------------------------

alg_latent <- function(F, S) {
  z <- if (S$nF == 2) F$t2("z\u{2c7c} = (z\u{2081}\u{2c7c}, z\u{2082}\u{2c7c})", "z_j = (z_{1j}, z_{2j})") else F$t2("z\u{2c7c}", "z_j")
  f1 <- S$fs[1]
  lines <- list(al_eq(F, z, F$sim, F$t2(sprintf("N(0, I%s)", if (S$nF == 2) "\u{2082}" else "\u{2081}"), sprintf("\\mathcal{N}(0, I_%d)", S$nF)),
                      "independent standard normals per person"),
                al_eq(F, F$s(F$v(f1), "j"), F$eq, F$t2("x\u{2c7c}\u{2032}\u{03b2}\u{2081} + s\u{2081}(x\u{2c7c}) z\u{2081}\u{2c7c}", "x_j'\\beta_1 + s_1(x_j)\\, z_{1j}"), "lat_values()"))
  if (S$nF == 2) {
    f2 <- S$f[[2]]
    Tu <- if (f2$shash) F$t2("T(u\u{2c7c})", "T(u_j)") else F$t2("u\u{2c7c}", "u_j")
    lines[[3]] <- al_eq(F, F$s(F$v(S$fs[2]), "j"), F$eq,
                        F$t2(sprintf("%sx\u{2c7c}\u{2032}\u{03b2}\u{2082} + s\u{2082} %s,  u\u{2c7c} = r z\u{2081}\u{2c7c} + \u{221a}(1 \u{2212} r\u{00b2}) z\u{2082}\u{2c7c}", if (f2$path) paste0("\u{03b3} ", F$s(F$v(f1), "j"), " + ") else "", Tu),
                             sprintf("%sx_j'\\beta_2 + s_2\\, %s,\\ u_j = r z_{1j} + \\sqrt{1 - r^2}\\, z_{2j}", if (f2$path) paste0("\\gamma ", F$s(F$v(f1), "j"), " + ") else "", Tu)),
                        if (f2$shash) "T: standardized sinh-arcsinh transform (shash())" else "")
  }
  lines[[length(lines) + 1]] <- al_txt(F, paste(F$t2("\u{03d1}:", "theta:"), "the working parameters (log of positive loadings and SDs, atanh of the correlation, ...); names(fit$rtmb$x)"))
  if (length(S$O)) lines[[length(lines) + 1]] <- al_txt(F, sprintf("ordinal indicators (%d): graded response, P(y > k) = logistic{a(f - b_k - delta)}, thresholds b_1 free, then b_k = b_(k-1) + exp(step); missing responses skipped",
                                                                   length(S$O)), "ind_ll()")
  if (!is.null(S$mod)) lines[[length(lines) + 1]] <- al_txt(F, sprintf("moderation: loadings times exp(alpha'z), shifts delta = beta'z (%d effects)", nrow(S$mod)), "mod_eff()")
  if (!is.null(S$cens)) lines[[length(lines) + 1]] <- al_txt(F, sprintf("censored indicators (%s; %d values): a value at or above its upper limit c contributes log S(c | f) = log Phi(-T^-1((c - mu) / sd)), at or below its lower limit log F(c | f) = log Phi(T^-1((c - mu) / sd)), instead of the log density (T: identity for normal, standardized SHASH transform)",
                                                                   paste(S$C[S$cens$j], collapse = ", "), as.integer(sum(S$Lc + S$Rc))), "ind_ll(), cens_logp()")
  list(title = "Person latents as transforms of standard normals", lines = lines)
}

alg_marginal <- function(F) list(title = "Marginal likelihood (the persons integrated out)", lines = list(
  al_eq(F, F$t2("L(\u{03d1})", "L(\\vartheta)"), F$eq,
        F$t2("\u{220f}\u{2c7c} \u{222b} p(y\u{2c7c} | f(z; \u{03d1}), \u{03d1}) \u{03c6}(z) dz", "\\prod_j \\int p(y_j \\mid f(z; \\vartheta), \\vartheta)\\, \\phi(z)\\, dz"),
        "one low-dimensional integral per person")))

alg_laplace <- function(F, S, staged, target = "maximum likelihood", prior = FALSE) {
  hz <- F$t2("\u{1e91}\u{2c7c}", "\\hat z_j"); Hi <- F$t2("H\u{2c7c}", "H_j")
  list(title = sprintf("%sLaplace approximation (%s)", if (target == "maximum likelihood") "" else "Start: ", target), lines = Filter(Negate(is.null), list(
    al_eq(F, F$t2("log p(y\u{2c7c} | \u{03d1})", "\\log p(y_j \\mid \\vartheta)"), F$t2("\u{2248}", "\\approx"),
          F$t2(sprintf("log p(y\u{2c7c}, %s | \u{03d1}) + (d/2) log 2\u{03c0} \u{2212} \u{00bd} log |%s|", hz, Hi),
               sprintf("\\log p(y_j, %s \\mid \\vartheta) + \\tfrac{d}{2}\\log 2\\pi - \\tfrac12 \\log|%s|", hz, Hi)),
          "RTMB::MakeADFun(random = z): inner Newton for the person modes"),
    al_eq(F, hz, F$eq, F$t2(sprintf("argmax_z log p(y\u{2c7c}, z | \u{03d1}),  %s = \u{2212}\u{2207}\u{00b2} log p(y\u{2c7c}, z | \u{03d1}) at %s", Hi, hz),
                            sprintf("\\arg\\max_z \\log p(y_j, z \\mid \\vartheta),\\ %s = -\\nabla^2 \\log p(y_j, z \\mid \\vartheta)\\big|_{%s}", Hi, hz)), ""),
    if (staged) al_txt(F, "stage 1: normal and homoscedastic (shape and moderation parameters fixed at 0); stage 2: the full model from there"),
    al_eq(F, F$t2("\u{03d1}\u{207d}\u{2070}\u{207e}", "\\vartheta^{(0)}"), F$eq,
          F$t2(sprintf("argmax \u{2211}\u{2c7c} log p\u{0303}(y\u{2c7c} | \u{03d1})%s", if (prior) " + log p(\u{03d1})" else ""),
               sprintf("\\arg\\max \\sum_j \\log \\tilde p(y_j \\mid \\vartheta)%s", if (prior) " + \\log p(\\vartheta)" else "")), "stats::nlminb with RTMB gradients"))))
}

alg_aghq <- function(F, S, k, rounds = NULL, ctrl = NULL, qchk = NULL, prior = FALSE) {
  d <- S$nF; G <- k^d
  zq <- F$t2("z\u{2c7c}\u{2099}", "z_{jn}")
  pr <- if (prior) F$t2(" + log p(\u{03d1})", " + \\log p(\\vartheta)") else ""
  centres <- list(title = "Quadrature centres per person", lines = list(
    al_eq(F, F$t2("(\u{1e91}\u{2c7c}, H\u{2c7c})", "(\\hat z_j, H_j)"), F$eq, F$t2("mode and negative Hessian of log p(y\u{2c7c}, z | \u{03d1}) in z", "\\text{mode and negative Hessian of } \\log p(y_j, z \\mid \\vartheta) \\text{ in } z"),
          "centers(): from the Laplace object at the current parameters"),
    al_eq(F, F$t2("L\u{2c7c} L\u{2c7c}\u{2032}", "L_j L_j'"), F$eq, F$t2("H\u{2c7c}\u{207b}\u{00b9}", "H_j^{-1}"), "Cholesky factor (2 x 2 in closed form)")))
  quad_l <- list(
    al_eq(F, zq, F$eq, F$t2("\u{1e91}\u{2c7c} + \u{221a}2 L\u{2c7c} x\u{2099}", "\\hat z_j + \\sqrt2\\, L_j x_n"),
          sprintf("x_n, w_n: %d Gauss-Hermite nodes per dimension (statmod::gauss.quad), %s", k, if (d == 2) sprintf("product rule, %d nodes per person", G) else sprintf("%d nodes per person", G))),
    al_eq(F, F$t2("log p(y\u{2c7c} | \u{03d1})", "\\log p(y_j \\mid \\vartheta)"), F$t2("\u{2248}", "\\approx"),
          F$t2(sprintf("log \u{2211}\u{2099} w\u{2099} e^{|x\u{2099}|\u{00b2}} 2^{d/2} |L\u{2c7c}| p(y\u{2c7c} | f(%s), \u{03d1}) \u{03c6}(%s)", zq, zq),
               sprintf("\\log \\sum_n w_n e^{\\lVert x_n \\rVert^2} 2^{d/2} |L_j|\\, p(y_j \\mid f(%s), \\vartheta)\\, \\phi(%s)", zq, zq)),
          "make_aghq(): the nodes are data, so this is an ordinary RTMB objective (log-sum-exp per person)"))
  quad <- list(title = sprintf("Adaptive Gauss-Hermite quadrature (k = %d)", k), lines = quad_l)
  if (is.null(ctrl)) return(list(centres, quad))
  nr <- length(rounds) - 1
  obj <- if (prior) "-2 log L - 2 log prior (the posterior mode objective)" else "-2 log L"
  opt <- list(title = if (prior) "Posterior mode: optimization rounds" else "Optimization rounds", lines = Filter(Negate(is.null), list(
    al_eq(F, F$t2("\u{03d1}\u{207d}\u{02b3}\u{207a}\u{00b9}\u{207e}", "\\vartheta^{(r+1)}"), F$eq,
          F$t2(sprintf("argmax log L_AGHQ(\u{03d1}; centres at \u{03d1}\u{207d}\u{02b3}\u{207e})%s", pr), sprintf("\\arg\\max \\log L_{\\mathrm{AGHQ}}(\\vartheta; \\text{centres at } \\vartheta^{(r)})%s", pr)),
          "stats::nlminb, RTMB gradients"),
    al_txt(F, "recentre the nodes at the new estimates and re-evaluate; if the objective got worse, take a half step (1/2, 1/4, 1/8, 1/16)"),
    al_txt(F, sprintf("stop when %s changes by less than tol = %s, at most max_rounds = %d; this fit: %d round%s",
                      obj, format(ctrl$tol), ctrl$max_rounds, nr, if (nr == 1) "" else "s")),
    al_txt(F, "polish at the final nodes, recentre once more and evaluate there"),
    if (!is.null(qchk)) al_txt(F, sprintf("accuracy check: %s = %s with k = %d and %s with k = %d (a warning if they differ by more than 0.5)",
                                          obj, format(round(qchk[1], 3), nsmall = 3), k, format(round(qchk[2], 3), nsmall = 3), k + 6)))))
  list(centres, quad, opt)
}

alg_se <- function(F, laplace = FALSE) list(title = "Standard errors", lines = list(
  if (laplace) al_eq(F, "V", F$eq, F$t2("(\u{2212}\u{2207}\u{00b2} log L\u{0303}(\u{03d1}\u{0302}))\u{207b}\u{00b9}", "\\left(-\\nabla^2 \\log \\tilde L(\\hat\\vartheta)\\right)^{-1}"), "RTMB::sdreport() of the Laplace objective")
  else al_eq(F, "V", F$eq, F$t2("(\u{2212}\u{2207}\u{00b2} log L_AGHQ(\u{03d1}\u{0302}))\u{207b}\u{00b9}", "\\left(-\\nabla^2 \\log L_{\\mathrm{AGHQ}}(\\hat\\vartheta)\\right)^{-1}"),
             "Hessian by automatic differentiation (obj$he); parameters at a bound are held fixed"),
  al_eq(F, F$t2("se(h(\u{03d1}\u{0302}))", "\\mathrm{se}(h(\\hat\\vartheta))"), F$eq, F$t2("\u{221a}(\u{2207}h\u{2032} V \u{2207}h)", "\\sqrt{\\nabla h' V \\nabla h}"),
        "delta method for every reported quantity (RTMB::ADREPORT Jacobian); also quantiles()")))

alg_scores <- function(F, mixture = FALSE) list(title = "Person scores", lines = list(
  al_eq(F, F$t2("E(f\u{2c7c} | y\u{2c7c})", "E(f_j \\mid y_j)"), F$t2("\u{2248}", "\\approx"),
        F$t2("\u{2211}\u{2099} \u{03c0}\u{2c7c}\u{2099} f(z\u{2c7c}\u{2099}),  \u{03c0}\u{2c7c}\u{2099} \u{221d} quadrature weight \u{00d7} p(y\u{2c7c}, z\u{2c7c}\u{2099} | \u{03d1})",
             "\\sum_n \\pi_{jn} f(z_{jn}),\\ \\pi_{jn} \\propto \\text{weight} \\times p(y_j, z_{jn} \\mid \\vartheta)"),
        if (mixture) "EAP and posterior SD per outer node, mixed with the outer node probabilities" else "EAP and posterior SD at the final nodes; reliability = var(EAP) / (var(EAP) + mean PSD^2)")))

alg_elgm <- function(F, object) {
  r <- object$rtmb; E <- r$elgm; S <- object$S
  g <- F$t2("\u{03d1}\u{2097}", "\\vartheta_l")
  nh <- if (E$nq) sprintf("%s (%d)", paste(E$hyper, collapse = ", "), E$nq) else "none"
  sc_lines <- if (E$nq && !is.null(E$scaling)) lapply(seq_len(nrow(E$scaling)), function(i)
    al_txt(F, sprintf("scaling of %s: %.2f below, %.2f above the mode", E$scaling$par[i], E$scaling$neg[i], E$scaling$pos[i])))
  list(
    list(title = "Priors on the working scale", lines = list(
      al_eq(F, F$t2("log p(\u{03d1})", "\\log p(\\vartheta)"), F$eq,
            F$t2("\u{2211}\u{2098} [log p\u{2098}(h\u{2098}(\u{03d1}\u{2098})) + log |h\u{2098}\u{2032}(\u{03d1}\u{2098})|],  h\u{2098} \u{2208} {identity, exp, tanh}",
                 "\\sum_m \\left[\\log p_m(h_m(\\vartheta_m)) + \\log|h_m'(\\vartheta_m)|\\right],\\ h_m \\in \\{\\mathrm{id}, \\exp, \\tanh\\}"),
            sprintf("rtmb_priors(): %d prior terms on the natural scale (equations(fit)), with the Jacobian", S$prior$n)))),
    list(title = "Hessian at the posterior mode", lines = list(
      al_eq(F, "V", F$eq, F$t2("(\u{2212}\u{2207}\u{00b2} log p(\u{03d1} | y) at \u{03d1}\u{0302})\u{207b}\u{00b9}", "\\left(-\\nabla^2 \\log p(\\vartheta \\mid y)\\big|_{\\hat\\vartheta}\\right)^{-1}"), "automatic differentiation"),
      al_txt(F, sprintf("split the parameters into quadrature directions %s and %d Gaussian directions", nh, E$n_gauss)))),
    list(title = sprintf("Outer quadrature over the parameters (%d nodes)", length(E$prob)), lines = c(list(
      al_eq(F, F$t2("L L\u{2032}", "L L'"), F$eq, "V", "Cholesky factor with the quadrature directions first (their spread is marginal)"),
      if (E$nq) al_eq(F, g, F$eq, F$t2("\u{03d1}\u{0302} + L (d \u{2218} z\u{2097}, 0)", "\\hat\\vartheta + L\\,(d \\circ z_l, 0)"),
                      sprintf("z_l: product grid of %d Gauss-Hermite nodes on %d direction%s (statmod::gauss.quad.prob)", E$k, E$nq, if (E$nq == 1) "" else "s"))
      else al_txt(F, "no quadrature directions: one node at the mode, Gaussian (Laplace) in every direction"),
      if (E$nq) al_eq(F, F$t2("d\u{2098}\u{00b1}", "d_m^{\\pm}"), F$eq,
                      F$t2("min(max(\u{221a}(2 / \u{0394}\u{2098}\u{00b1}), 1/3), 3),  \u{0394}\u{2098}\u{00b1} = drop of log p(\u{03d1} | y) at z\u{2098} = \u{00b1}2",
                           "\\min(\\max(\\sqrt{2 / \\Delta_m^{\\pm}}, 1/3), 3),\\ \\Delta_m^{\\pm} = \\text{drop of } \\log p(\\vartheta \\mid y) \\text{ at } z_m = \\pm 2"),
                      "asymmetric scaling per side, as INLA's grid (1 for a Gaussian posterior)"),
      if (E$cond_modes) al_txt(F, "Gaussian directions at each node: their conditional posterior mode (Newton steps with the conditional covariance, step halving)", "cond_modes = TRUE")
      else if (E$nq && E$n_gauss) al_txt(F, "Gaussian directions at each node: on the Gaussian conditional-mean line of the mode", "cond_modes = FALSE")),
      sc_lines)),
    list(title = "Node weights and the marginal likelihood", lines = list(
      al_eq(F, F$t2("log p(y, \u{03d1}\u{2097})", "\\log p(y, \\vartheta_l)"), F$eq,
            F$t2("\u{2211}\u{2c7c} log p(y\u{2c7c} | \u{03d1}\u{2097}) + log p(\u{03d1}\u{2097})", "\\sum_j \\log p(y_j \\mid \\vartheta_l) + \\log p(\\vartheta_l)"),
            sprintf("inner %s per person at each node (aghq_value())", if (E$inner == "aghq") sprintf("AGHQ, k = %d", E$k_inner) else "Laplace")),
      al_eq(F, F$t2("log w\u{2097}", "\\log w_l"), F$eq,
            F$t2("log \u{03c9}\u{2097} \u{2212} log \u{03c6}(z\u{2097}) + log J\u{2097} + log |L| + (n_G / 2) log 2\u{03c0}",
                 "\\log \\omega_l - \\log \\phi(z_l) + \\log J_l + \\log|L| + \\tfrac{n_G}{2}\\log 2\\pi"),
            sprintf("J_l: Jacobian of the scaling; n_G = %d Gaussian directions", E$n_gauss)),
      al_eq(F, F$t2("log p(y)", "\\log p(y)"), F$t2("\u{2248}", "\\approx"),
            F$t2(sprintf("log \u{2211}\u{2097} w\u{2097} p(y, \u{03d1}\u{2097}) = %.2f", E$lognc), sprintf("\\log \\sum_l w_l\\, p(y, \\vartheta_l) = %.2f", E$lognc)),
            "fit_indices(): logML (Bayes factors in compare())"),
      al_eq(F, F$t2("\u{03c0}\u{2097}", "\\pi_l"), F$eq, F$t2("w\u{2097} p(y, \u{03d1}\u{2097}) / p(y)", "w_l\\, p(y, \\vartheta_l) / p(y)"),
            sprintf("posterior node probabilities; ESS of the nodes %.1f", E$ess)))),
    list(title = "Posterior summaries", lines = list(
      al_eq(F, F$t2("E(h | y)", "E(h \\mid y)"), F$t2("\u{2248}", "\\approx"), F$t2("\u{2211}\u{2097} \u{03c0}\u{2097} h(\u{03d1}\u{2097})", "\\sum_l \\pi_l\\, h(\\vartheta_l)"), "every reported quantity h (estimates(), quantiles())"),
      al_eq(F, F$t2("Var(h | y)", "\\mathrm{Var}(h \\mid y)"), F$t2("\u{2248}", "\\approx"),
            F$t2("\u{2211}\u{2097} \u{03c0}\u{2097} (h(\u{03d1}\u{2097}) \u{2212} E h)\u{00b2} + \u{2207}h\u{2032} \u{03a3}_c \u{2207}h", "\\sum_l \\pi_l (h(\\vartheta_l) - E h)^2 + \\nabla h' \\Sigma_c \\nabla h"),
            "Sigma_c: conditional covariance of the Gaussian directions given the quadrature ones"),
      al_txt(F, "95% interval: Cornish-Fisher from the quadrature mean, SD and skewness, on the log (SDs, positive loadings), atanh (correlations) or identity scale", "elgm_summary()"))))
}

#' @rdname algorithm
#' @export
algorithm.birt_rtmb <- function(object, format = c("unicode", "latex"), ...) {
  format <- match.arg(format); F <- alg_fmt(format)
  S <- object$S; r <- object$rtmb; ctrl <- object$control
  staged <- length(intersect(fam_pars, names(Filter(length, S$par0)))) > 0
  if (!is.null(r$elgm)) {
    E <- r$elgm; inner_aghq <- r$method == "aghq"
    secs <- c(list(alg_latent(F, S), alg_marginal(F)), alg_elgm(F, object)[1],
              list(alg_laplace(F, S, staged, "posterior mode", prior = TRUE)),
              if (inner_aghq) alg_aghq(F, S, r$k, r$rounds, ctrl, r$quad_check, prior = TRUE)
              else list(list(title = "Inner integral: Laplace", lines = list(al_txt(F, "the Laplace value above is used for log p(y_j | theta) at the mode and at every outer node", "rtmb_control(inner = \"laplace\")")))),
              alg_elgm(F, object)[-1], list(alg_scores(F, mixture = TRUE)))
    return(alg_object(secs, format))
  }
  if (r$method == "laplace") {
    secs <- c(list(alg_latent(F, S), alg_marginal(F), alg_laplace(F, S, staged)), list(alg_se(F, laplace = TRUE)),
              list(list(title = sprintf("Person scores (AGHQ with k = %d at the Laplace centres)", r$k), lines = alg_scores(F)$lines)))
    return(alg_object(secs, format))
  }
  secs <- c(list(alg_latent(F, S), alg_marginal(F), alg_laplace(F, S, staged, "start values")),
            alg_aghq(F, S, r$k, r$rounds, ctrl, r$quad_check), list(alg_se(F), alg_scores(F)))
  alg_object(secs, format)
}

