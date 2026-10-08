# Model specification: syntax + data checks + family ------------------------------

build_spec <- function(model, data, binary = NULL, family = NULL, ordered = NULL, itemtype = "grm") {
  if (!is.data.frame(data)) stop("`data` must be a data frame")
  sp <- parse_model(model)
  vars <- sp$ovs
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("variables in the model but not in `data`: ", paste(miss, collapse = ", "))
  nonnum <- vars[!vapply(vars, function(v) is.numeric(data[[v]]) || is.logical(data[[v]]), logical(1))]
  if (length(nonnum)) stop("non-numeric indicators: ", paste(nonnum, collapse = ", "))
  allna <- vars[vapply(vars, function(v) all(is.na(data[[v]])), logical(1))]
  if (length(allna)) stop("indicators with no observed values: ", paste(allna, collapse = ", "))
  const <- vars[vapply(vars, function(v) length(unique(stats::na.omit(data[[v]]))) == 1, logical(1))]
  if (length(const))
    stop(sprintf("indicators without variation (one observed value): %s; their parameters are not identified (e.g. an item everyone answers correctly), drop them from the model",
                 paste(sprintf("%s = %s", const, vapply(const, function(v) format(stats::na.omit(data[[v]])[1]), "")), collapse = ", ")), call. = FALSE)
  # ordinal indicators: integer categories, recoded to 1..K per indicator
  isint <- vapply(vars, function(v) { x <- stats::na.omit(data[[v]]); is.numeric(x) && all(x == round(x)) }, TRUE)
  ncat <- vapply(vars, function(v) length(unique(stats::na.omit(data[[v]]))), 0L)
  ord <- if (is.null(ordered) || isFALSE(ordered)) character() else if (isTRUE(ordered)) {
    o <- vars[isint & ncat >= 3 & ncat <= 10]
    other <- vars[isint & !vars %in% o & ncat > 2]
    if (length(other)) message("ordered = TRUE: integer indicators with more than 10 values stay continuous (", paste(utils::head(other, 4), collapse = ", "), "); name them in `ordered` to make them ordinal")
    o } else ordered
  if (!is.null(itemtype) && !identical(unname(itemtype), "grm") && !length(ord))
    stop("itemtype applies to ordinal indicators: give `ordered`", call. = FALSE)
  if (!is.character(ord)) stop("`ordered` must be TRUE or the names of the ordinal indicators", call. = FALSE)
  bad <- setdiff(ord, vars); if (length(bad)) stop("`ordered` variables not in the model: ", paste(bad, collapse = ", "), call. = FALSE)
  bad <- ord[!isint[ord]]; if (length(bad)) stop("ordinal indicators need integer categories: ", paste(bad, collapse = ", "), call. = FALSE)
  sp$ord <- list(items = character(), levels = list(), K = integer(), Y = NULL, type = character())   # exists always: no partial matching of sp$ord
  if (length(ord)) {
    lv <- lapply(stats::setNames(ord, ord), function(v) sort(unique(stats::na.omit(data[[v]]))))
    Y <- vapply(ord, function(v) match(data[[v]], lv[[v]]), integer(nrow(data))); Y <- matrix(Y, nrow(data), dimnames = list(NULL, ord))
    K <- vapply(lv, length, 1L)
    short <- names(K)[K < max(K)]
    if (length(short)) message("items with fewer observed categories (recoded to 1..K): ",
                               paste(sprintf("%s (%s)", short, vapply(lv[short], paste, "", collapse = "/")), collapse = ", "))
    if (!is.character(itemtype) || !all(itemtype %in% ord_types)) stop("itemtype must be one of: ", paste(ord_types, collapse = ", "), call. = FALSE)
    ty <- if (!is.null(names(itemtype))) { bad <- setdiff(names(itemtype), ord)
      if (length(bad)) stop("itemtype: not ordinal indicators of the model: ", paste(bad, collapse = ", "), call. = FALSE)
      x <- stats::setNames(rep("grm", length(ord)), ord); x[names(itemtype)] <- itemtype; x }
          else if (length(itemtype) == 1) stats::setNames(rep(itemtype, length(ord)), ord)
          else stop("itemtype: one value, or a named vector (item = type)", call. = FALSE)
    sp$ord <- list(items = ord, levels = lv, K = unname(K), Y = Y, type = ty)
  }
  is01 <- vapply(vars, function(v) all(data[[v]] %in% c(0, 1, NA)), logical(1)) & !vars %in% ord
  two <- vars[!is01 & !vars %in% ord & vapply(vars, function(v) length(unique(stats::na.omit(data[[v]]))) == 2, logical(1))]
  if (length(two)) {
    lv <- sort(unique(stats::na.omit(data[[two[1]]])))
    stop(sprintf("two-valued indicators not coded 0/1 (%s: %s/%s%s); binary items must be 0/1, e.g. data$%s <- data$%s - %s",
                 two[1], lv[1], lv[2], if (length(two) > 1) sprintf("; also %s", paste(two[-1], collapse = ", ")) else "",
                 two[1], two[1], lv[1]), call. = FALSE)
  }
  if (is.null(binary)) {
    binary <- vars[is01]
  } else {
    if (length(intersect(binary, ord))) stop("indicators both binary and ordered: ", paste(intersect(binary, ord), collapse = ", "), call. = FALSE)
    bad <- setdiff(binary, vars)
    if (length(bad)) stop("`binary` variables not in the model: ", paste(bad, collapse = ", "))
    notbin <- binary[!is01[binary]]
    if (length(notbin)) stop("binary indicators with values other than 0/1: ", paste(notbin, collapse = ", "))
  }
  sp$binary <- vars[vars %in% binary]
  sp$ordinal <- vars[vars %in% ord]
  sp$cont <- setdiff(vars, c(binary, ord))
  raw_rt <- sp$cont[vapply(sp$cont, function(v) {
    x <- stats::na.omit(data[[v]])
    length(x) > 2 && all(x > 0) && stats::median(x) > 10 && mean((x - mean(x))^3) / stats::sd(x)^3 > 1
  }, logical(1))]
  if (length(raw_rt))
    warning(sprintf("continuous indicators look like raw response times (positive, right-skewed, median > 10: %s%s); the model assumes normal residuals, so log-transform them: birt(..., log_rt = TRUE)%s%s",
                    paste(utils::head(raw_rt, 4), collapse = ", "), if (length(raw_rt) > 4) ", ..." else "",
                    "", ""),
            call. = FALSE)
  ratio <- vapply(vars, function(v) mean(is.na(data[[v]])), 0)
  if (any(ratio > 0.5)) warning("more than half missing: ", paste(vars[ratio > 0.5], collapse = ", "))

  miss_cov <- setdiff(sp$covs, names(data))
  if (length(miss_cov)) stop("covariates in the model but not in `data`: ", paste(miss_cov, collapse = ", "))
  na_cov <- sp$covs[vapply(sp$covs, function(v) anyNA(data[[v]]) || !is.numeric(data[[v]]), logical(1))]
  if (length(na_cov)) stop("covariates must be numeric without missing values: ", paste(na_cov, collapse = ", "))

  # identification: a correlation / path between two latent variables together with free
  # cross-loadings of one on every indicator (with fixed loadings) of the other
  sp$ident <- NULL
  if (length(sp$factors) == 2) {
    ld <- sp$ld[!sp$ld$zero, ]
    for (o in list(sp$factors, rev(sp$factors))) {
      a <- o[1]; b <- o[2]
      ib <- ld$rhs[ld$lhs == b]
      link <- isTRUE(sp$cov_free) || any(sp$reg$is_factor & ((sp$reg$lhs == b & sp$reg$x == a) | (sp$reg$lhs == a & sp$reg$x == b)))
      cr <- ld[ld$lhs == a & ld$rhs %in% ib & ld$fixed == "", ]
      if (link && length(ib) && all(ld$fixed[ld$lhs == b] != "") && all(ib %in% cr$rhs))
        sp$ident <- list(a = a, b = b, prior = all(cr$prior[match(ib, cr$rhs)] != ""),
                         msg = sprintf("%s has cross-loadings on every indicator of %s, whose loadings are fixed, together with %s: a common cross-loading is the same as the %s, so they are not identified by the likelihood; drop the %s or anchor some cross-loadings (e.g. 0*%s)",
                                       a, b, if (isTRUE(sp$cov_free)) sprintf("%s ~~ %s", sp$factors[1], sp$factors[2]) else "the latent regression path",
                                       if (isTRUE(sp$cov_free)) "correlation" else "path", if (isTRUE(sp$cov_free)) "correlation" else "path", ib[1]))
    }
  }
  sp$family <- resolve_family(family, sp)
  sp <- check_moderation(sp, data)
  sp
}

# The moderation statements against the model and the data. V() statements become the scale
# moderators of the families (the log SD / log scale of the residual distribution).
check_moderation <- function(sp, data) {
  mo <- sp$mod; ld <- sp$ld[!sp$ld$zero, ]; fs <- sp$factors
  cat_items <- c(sp$binary, sp$ordinal)
  for (v in sp$ordinal) {
    r <- ld[ld$rhs == v, ]
    if (nrow(r) != 1) stop("ordinal indicator ", v, " must load on exactly one latent variable", call. = FALSE)
    if (if (r$fixed != "") !(as.numeric(r$fixed) > 0) else !r$first)
      stop("ordinal indicator ", v, ": its loading (discrimination) must be positive: a free loading on the first line of ", r$lhs, " or a fixed positive value", call. = FALSE)
  }
  obs <- setdiff(unique(mo$mod), fs)
  miss <- setdiff(obs, names(data)); if (length(miss)) stop("moderators not in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  bad <- obs[vapply(obs, function(v) !is.numeric(data[[v]]) || anyNA(data[[v]]), TRUE)]
  if (length(bad)) stop("moderators must be numeric without missing values: ", paste(bad, collapse = ", "), call. = FALSE)
  for (i in seq_len(nrow(mo))) {
    r <- mo[i, ]; stmt <- sprintf("%s(%s) ~ %s%s", if (r$kind %in% c("psi", "kappa")) "V" else "E", r$target, r$mod, if (nzchar(r$factor)) paste0(":", r$factor) else "")
    if (r$kind %in% c("alpha", "beta") && r$mod %in% fs) stop(stmt, ": the moderator of a loading or a shift is an observed variable", call. = FALSE)
    if (r$kind == "alpha" && !any(ld$lhs == r$factor & ld$rhs == r$target)) stop(stmt, ": ", r$target, " does not load on ", r$factor, call. = FALSE)
    if (r$kind == "beta" && r$target %in% cat_items && sum(ld$rhs == r$target) != 1)
      stop(stmt, ": a shift of a categorical indicator is on the latent scale, so the indicator must load on one latent variable", call. = FALSE)
    if (r$kind == "kappa" && !r$target %in% sp$cont) stop(stmt, ": categorical indicators have no residual SD to moderate (logistic, fixed scale)", call. = FALSE)
  }
  # V() -> scale moderators of the family (with their priors)
  for (i in which(mo$kind %in% c("psi", "kappa"))) {
    r <- mo[i, ]; slot <- if (r$kind == "psi") "factor" else "ind"
    if (r$type == "none") next                                           # V(y) ~ 0*z: no effect
    if (r$type == "common") stop(sprintf("V(%s) ~ %s*%s: equal (labelled) effects in V() are not supported", r$target, r$label, r$mod), call. = FALSE)
    f <- sp$family[[slot]][[r$target]]
    f$scale <- c(f$scale, r$mod); f$scale_prior <- c(f$scale_prior %||% character(), stats::setNames(r$prior, r$mod))
    sp$family[[slot]][[r$target]] <- f
  }
  sp$mod <- mo[mo$kind %in% c("alpha", "beta"), ]
  # identification: a shift of every indicator of f on z together with E(f) ~ z
  for (f in fs) for (z in sp$reg$x[sp$reg$lhs == f]) {
    ind <- unique(ld$rhs[ld$lhs == f])
    sh <- sp$mod$target[sp$mod$kind == "beta" & sp$mod$mod == z & sp$mod$type %in% c("free", "common")]
    if (length(ind) && all(ind %in% sh))
      stop(sprintf("E(%s) ~ %s with shifts of every indicator of %s on %s (E(item) ~ %s) is not identified (a common shift equals a mean shift): leave at least one indicator out (an anchor), or use prior(\"ssp\")",
                   f, z, f, z, z), call. = FALSE)
  }
  for (f in fs) for (z in sp$family$factor[[f]]$scale) {
    ind <- unique(ld$rhs[ld$lhs == f])
    al <- sp$mod$target[sp$mod$kind == "alpha" & sp$mod$factor == f & sp$mod$mod == z & sp$mod$type %in% c("free", "common")]
    if (length(ind) && all(ind %in% al))
      message("V(", f, ") ~ ", z, " with loading effects of every indicator: the common part of the loading effect and the variance effect are separated only through the intercepts / thresholds (weakly identified)")
  }
  sp$modvars <- obs
  sp
}
