# Model syntax ---------------------------------------------------------------------
#
# lavaan operators =~, ~, ~~ with modifiers (fixed values, labels, prior("...")), parsed
# by lavaan's parser, and the moderation statements of moderated nonlinear factor analysis
# (Bauer, 2017), which regress a conditional moment on moderators z:
#   E(y) ~ z      shift of indicator y: categorical y on the latent scale, a(f - b - delta),
#                 delta = beta'z (uniform DIF); continuous y its intercept
#   E(y) ~ z:f    loading of y on f, log scale: lambda(z) = lambda exp(alpha'z) (nonuniform DIF)
#   V(y) ~ z      log residual SD of a continuous y
#   E(f) ~ z      latent mean, the same as the regression f ~ z
#   V(f) ~ z      log SD of the latent variable f

parse_model <- function(model, binary = NULL) {
  if (!is.character(model) || length(model) != 1) stop("`model` must be a single character string")
  model <- as.character(unclass(model))
  check_modifiers(model)
  pt <- lavaan::lavParseModelString(mod_preparse(model), as.data.frame. = TRUE, parser = "old")
  check_labels(pt)
  cons <- attr(pt, "constraints")
  if (length(cons)) stop("defined parameters and constraints (", paste(unique(vapply(cons, `[[`, "", "op")), collapse = " "),
                         ") are not supported; use labels for equality constraints", call. = FALSE)
  bad_op <- setdiff(unique(pt$op), c("=~", "~", "~~"))
  if (length(bad_op)) stop("unsupported operator(s): ", paste(bad_op, collapse = " "))
  ld <- pt[pt$op == "=~", ]
  if (!nrow(ld)) stop("the model needs a measurement part (=~)")
  factors <- unique(ld$lhs)
  if (length(factors) > 2) stop("birt supports at most two latent variables (found ", length(factors), ")")
  ovs <- unique(ld$rhs)
  if (any(ovs %in% factors)) stop("higher-order factors are not supported")
  ld$line <- cumsum(c(TRUE, ld$lhs[-1] != ld$lhs[-nrow(ld)]))
  ld$first <- stats::ave(ld$line, ld$lhs, FUN = function(x) x == min(x)) == 1
  ld$zero <- ld$fixed != "" & suppressWarnings(as.numeric(ld$fixed)) %in% 0
  if (any(ld$prior != "" & ld$fixed != "")) stop("a loading cannot have both a fixed value and a prior")

  rg <- pt[pt$op == "~", ]
  mo <- rg[grepl("^birt[EV]__", rg$lhs), ]; rg <- rg[!grepl("^birt[EV]__", rg$lhs), ]
  mt <- mod_table(mo, factors, ovs)
  shared <- intersect(mt$mod$label[nzchar(mt$mod$label)], c(ld$label, rg$label)[nzchar(c(ld$label, rg$label))])
  if (length(shared)) stop(sprintf("label '%s' is shared by a moderation effect and a loading or regression coefficient; a label makes parameters of one kind equal", shared[1]), call. = FALSE)
  rg <- rbind(rg, mt$reg)                                              # E(f) ~ z is the regression f ~ z
  bad <- setdiff(rg$lhs, factors)
  if (length(bad)) stop("regression outcomes must be latent variables: ", paste(bad, collapse = ", "))
  rg$x <- rg$rhs
  if (any(rg$fixed != "")) stop("fixed regression coefficients are not supported")
  if (any(duplicated(paste(rg$lhs, rg$x))))
    stop("each predictor may enter a regression once")
  if (any(rg$x %in% ovs)) stop("indicators cannot be predictors in a latent regression")
  fpath <- rg[rg$x %in% factors, ]
  if (nrow(fpath)) {
    if (nrow(fpath) > 1 || fpath$lhs == fpath$x) stop("only one path between the two latent variables is supported")
    factors <- c(fpath$x, fpath$lhs)                     # predictor first
  }
  rg$is_factor <- rg$x %in% factors
  covs <- setdiff(unique(rg$x), factors)

  cv <- pt[pt$op == "~~" & pt$lhs != pt$rhs, ]
  if (nrow(cv) && !all(cv$lhs %in% factors & cv$rhs %in% factors))
    stop("residual covariances between indicators are not supported")
  if (nrow(cv) && any(cv$fixed != "" & !(suppressWarnings(as.numeric(cv$fixed)) %in% 0)))
    stop("a fixed nonzero factor covariance is not supported")
  cov_free <- nrow(cv) > 0 && any(cv$fixed == "")
  if (cov_free && nrow(fpath)) stop("use either 'f1 ~~ f2' or the regression 'f2 ~ f1', not both")
  iv <- pt[pt$op == "~~" & pt$lhs == pt$rhs & !(pt$lhs %in% factors), ]
  if (nrow(iv)) stop("indicator variances are always free; remove: ", paste(iv$lhs, "~~", iv$rhs, collapse = ", "))

  # factor (residual) SD: "f ~~ c*f" fixed, "f ~~ f" free, otherwise 1 when the factor
  # has a free loading and free when all its loadings are fixed
  fv <- pt[pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% factors, ]
  fsd <- vapply(factors, function(f) {
    r <- fv[fv$lhs == f, ]
    if (nrow(r) && r$fixed[1] != "") {
      v <- as.numeric(r$fixed[1]); if (!(v > 0)) stop("a fixed factor variance must be positive")
      return(sqrt(v))
    }
    if (nrow(r) || all(ld$fixed[ld$lhs == f] != "")) return(NA_real_)
    1
  }, numeric(1))

  # a path f2 ~ f1 is aliased with free loadings of f1 on all the indicators of f2 (with
  # fixed loadings on f2, a common cross-loading equals the path coefficient); cross-loadings
  # on some of them are identified (the others anchor the path)
  if (nrow(fpath)) {
    f1 <- factors[1]; f2 <- factors[2]
    ind2 <- ld$rhs[ld$lhs == f2 & !ld$zero]
    cross <- ld$rhs[ld$lhs == f1 & ld$rhs %in% ind2 & ld$fixed == ""]
    if (length(cross) && length(cross) == length(ind2))       # every indicator of f2: aliased
      stop("the path '", f2, " ~ ", f1, "' is not identified together with the free loadings of '", f1,
           "' on indicators of '", f2, "' (", paste(utils::head(cross, 3), collapse = ", "),
           if (length(cross) > 3) ", ..." else "", "); drop one of them")
  }
  list(pt = pt, ld = ld, factors = factors, ovs = ovs, reg = rg, covs = covs, fsd = fsd,
       cov_free = cov_free, binary = binary, cont = setdiff(ovs, binary), mod = mt$mod)
}

# moderation statements -> placeholder outcomes that lavaan's parser reads (birtE__y1 ~ z)
mod_preparse <- function(model) {
  lines <- unlist(strsplit(model, "[\n;]"))
  out <- vapply(lines, function(l) {
    l0 <- sub("#.*", "", l)
    m <- regmatches(l0, regexec("^\\s*([EV])\\s*\\(([^)]*)\\)\\s*~(.*)$", l0))[[1]]
    if (!length(m)) return(l0)
    v <- trimws(strsplit(m[3], "+", fixed = TRUE)[[1]])
    sprintf("%s ~%s", paste0("birt", m[2], "__", v, collapse = " + "), m[4])
  }, "")
  paste(out, collapse = "\n")
}

# The moderation statements as a table: kind (alpha: log loading, beta: shift, kappa: log residual
# SD, psi: log SD of a latent variable), target (indicator or latent variable), factor (alpha: the
# latent variable of the loading), mod (moderator), type (free, ssp, common = labelled, none = 0*z
# anchor), prior and label. E(f) ~ z rows come back as regressions.
# the empty moderation table (one row per moderation effect; see mod_table())
mod_empty <- function() data.frame(kind = character(), target = character(), factor = character(), mod = character(),
                                   type = character(), prior = character(), label = character(), stringsAsFactors = FALSE)

mod_table <- function(mo, factors, ovs) {
  empty <- mod_empty()
  if (!nrow(mo)) return(list(mod = empty, reg = mo))
  rows <- list(); reg <- mo[0, ]
  for (i in seq_len(nrow(mo))) {
    tag <- substr(mo$lhs[i], 5, 5); tg <- sub("^birt[EV]__", "", mo$lhs[i])
    parts <- strsplit(mo$rhs[i], ":", fixed = TRUE)[[1]]
    fpart <- parts[parts %in% factors]; z <- parts[!parts %in% factors]
    stmt <- sprintf("%s(%s) ~ %s", tag, tg, mo$rhs[i])
    if (length(parts) > 2 || length(z) > 1 || any(z %in% ovs))
      stop(stmt, ": a moderator is an observed variable (not an indicator), possibly times one latent variable (z:f)", call. = FALSE)
    if (!length(z)) {                                                  # V(f2) ~ f1: a latent moderator
      if (length(fpart) != 1 || tag != "V") stop(stmt, ": not a moderation statement", call. = FALSE)
      z <- fpart; fpart <- character()
    }
    inter <- length(fpart) == 1
    if (!(tg %in% c(factors, ovs))) stop(stmt, ": ", tg, " is not a latent variable or an indicator of the model", call. = FALSE)
    if (tg %in% factors && inter) stop(stmt, ": the mean and SD of a latent variable take moderators, not interactions", call. = FALSE)
    if (tag == "E" && tg %in% factors) {
      if (mo$fixed[i] != "") stop(stmt, ": fixed regression coefficients are not supported", call. = FALSE)
      r <- mo[i, ]; r$lhs <- tg; r$rhs <- z; reg <- rbind(reg, r); next
    }
    if (tag == "V" && inter) stop(stmt, ": the SD takes moderators, not interactions", call. = FALSE)
    if (mo$fixed[i] != "" && !isTRUE(suppressWarnings(as.numeric(mo$fixed[i])) == 0))
      stop(stmt, ": fixed nonzero moderation effects are not supported (0*z marks an anchor)", call. = FALSE)
    kind <- if (tag == "V") (if (tg %in% factors) "psi" else "kappa") else if (inter) "alpha" else "beta"
    type <- if (mo$fixed[i] != "") "none" else if (mo$prior[i] == "ssp") "ssp" else if (mo$label[i] != "") "common" else "free"
    if (type == "ssp" && mo$label[i] != "") stop(stmt, ": a labelled (equal) effect cannot have prior(\"ssp\")", call. = FALSE)
    if (type == "ssp" && kind %in% c("psi", "kappa")) stop(stmt, ": prior(\"ssp\") is for E() effects", call. = FALSE)
    rows[[length(rows) + 1]] <- data.frame(kind = kind, target = tg, factor = if (inter) fpart else "", mod = z, type = type,
                                           prior = if (type == "ssp") "" else mo$prior[i], label = mo$label[i], stringsAsFactors = FALSE)
  }
  mod <- if (length(rows)) do.call(rbind, rows) else empty
  if (any(duplicated(mod[c("kind", "target", "factor", "mod")]))) stop("a moderation term is given twice", call. = FALSE)
  for (g in unique(mod$label[mod$type == "common"])) {
    r <- mod[mod$type == "common" & mod$label == g, ]
    if (length(unique(r$kind)) > 1) stop("label ", g, " is shared by different kinds of moderation effects", call. = FALSE)
    if (length(unique(r$mod)) > 1) stop("label ", g, ": equal effects must belong to the same moderator", call. = FALSE)
    if (length(unique(r$prior[r$prior != ""])) > 1) stop("label ", g, ": different priors on equal effects", call. = FALSE)
  }
  list(mod = mod, reg = reg)
}

# One modifier per term: lavaan's parser keeps only the label of "0.5*p*y" or "(p*0.5)*y" and
# drops the number, so a scaled label would silently become plain equality.
check_modifiers <- function(model) {
  x <- gsub('"[^"]*"|\'[^\']*\'', '""', model)                        # quoted prior strings
  x <- gsub("#[^\n]*", "", x)                                          # comments
  x <- gsub("\\+\\s*\n\\s*", "+ ", gsub("\\s*\n\\s*\\+", " +", x))          # continuation lines
  x <- gsub("\\bprior\\([^()]*\\)", "P", x, perl = TRUE)                                    # prior("...")
  repeat { y <- gsub("\\b[A-Za-z_.][A-Za-z0-9_.]*\\([^()]*\\)", "F", x, perl = TRUE); if (y == x) break; x <- y }   # start(), c(), ...
  for (ln in unlist(strsplit(x, "[\n;]"))) {
    rhs <- sub("^.*?(=~|~~|~)", "", ln, perl = TRUE)
    if (identical(rhs, ln)) next
    for (tm in strsplit(rhs, "+", fixed = TRUE)[[1]]) {
      tm <- trimws(tm); if (!nzchar(tm)) next
      parts <- trimws(strsplit(tm, "*", fixed = TRUE)[[1]])
      mods <- parts[-length(parts)]
      lab <- mods[!mods %in% c("F", "P") & is.na(suppressWarnings(as.numeric(mods)))]
      if (grepl("[()]", tm) || sum(!mods %in% c("F", "P")) > 1)
        stop(sprintf("'%s': a term takes one value or one label (plus prior() / start()); a label times a number (e.g. (p*0.025)*t) is not supported", tm), call. = FALSE)
      if (length(lab) && any(mods == "P"))
        stop(sprintf("'%s': a label and prior() on one term: lavaan's parser keeps only the label (the prior would be lost silently); labelled parameters take the default prior (dpriors())", tm), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# A label shared by several parameters makes them equal. Sharing is supported among loadings
# and among regression coefficients; across the two, or on variances / covariances, it would
# be ignored, so it is refused.
check_labels <- function(pt) {
  lb <- pt[pt$label != "", ]
  if (!nrow(lb)) return(invisible(TRUE))
  shared <- names(which(table(lb$label) > 1))
  for (l in shared) {
    ops <- unique(lb$op[lb$label == l])
    if ("~~" %in% ops) stop(sprintf("label '%s' is shared by a variance or covariance; equality constraints on (co)variances are not supported", l), call. = FALSE)
    if (length(ops) > 1) stop(sprintf("label '%s' is shared by a loading (=~) and a regression coefficient (~); a label can make loadings equal, or regression coefficients equal, not both", l), call. = FALSE)
  }
  invisible(TRUE)
}
