# Summary objects: tables computed once, printed in the layout of tppcm::irt_pars() and
# blavaan --------------------------------------------------------------------------------
#
# summary(fit) returns a "birt_summary": a list of titled tables (data frames of class
# "birt_table") plus header lines. Every table can be extracted (s$items, s$parameters,
# as.data.frame(s)); print() shows them under their names ($items, $se, ...), the parameter
# table in lavaan's sections (Latent Variables, Regressions, ...) with a Prior column for
# Bayesian fits.

btable <- function(df, title, cls = NULL) {
  if (is.null(df)) return(NULL)
  rownames(df) <- NULL
  structure(df, class = c(cls, "birt_table", "data.frame"), title = title)
}

#' Summary tables
#'
#' `summary()` of a birt fit returns the tables of the fit as one object: `header`, `fit`
#' (fit indices), `items` and `se` (one row per item: estimates and their posterior SD or
#' standard error, as `coef(simplify = TRUE)` in mirt), `parameters` (one row per parameter
#' with its interval and, for Bayesian fits, its prior, as in blavaan) and model-specific
#' tables (e.g. `selection`, `moderation`, `precision`, `reliability`). Each element is a data
#' frame; `as.data.frame()` returns `parameters`. `print()` shows the tables under their
#' names, so `s$items` gets the one printed as `$items`.
#' @param x A `birt_summary`, or one of its tables.
#' @param tables Names of the tables to print (default all).
#' @param digits Decimals.
#' @param title Print the title of a table.
#' @param row.names,optional Unused.
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `as.data.frame()` the parameter table.
#' @name birt_summary
NULL

#' @rdname birt_summary
#' @export
print.birt_summary <- function(x, tables = NULL, digits = 3, ...) {
  cat(x$header, sep = "\n")
  nm <- setdiff(names(x), "header")
  for (n in if (is.null(tables)) nm else intersect(tables, nm)) {
    t <- x[[n]]
    if (is.null(t) || (is.data.frame(t) && !nrow(t))) next
    cat(sprintf("\n$%s%s\n", n, if (!is.null(attr(t, "title"))) paste0("  ", attr(t, "title")) else ""))
    print(t, digits = digits, title = FALSE)
  }
  invisible(x)
}

#' @rdname birt_summary
#' @export
as.data.frame.birt_summary <- function(x, row.names = NULL, optional = FALSE, ...) {
  p <- x$parameters; attr(p, "title") <- NULL; class(p) <- "data.frame"; p
}

#' @rdname birt_summary
#' @export
print.birt_table <- function(x, digits = 3, title = TRUE, ...) {
  if (title && !is.null(attr(x, "title"))) cat(attr(x, "title"), "\n")
  d <- x; attr(d, "title") <- NULL; class(d) <- "data.frame"
  print_table(d, digits)
  invisible(x)
}

# ---- the parameter table: lavaan sections, blavaan columns ------------------------------

# standard columns from an estimates() table: est, se (posterior SD or standard error),
# lower, upper, then rhat, ess, p_incl, BF10, prior (Bayes) or z, pvalue (ML)
std_partable <- function(e, bayes) {
  g <- function(n) if (n %in% names(e)) e[[n]] else rep(NA, nrow(e))
  out <- data.frame(lhs = e$lhs, op = e$op, rhs = e$rhs, est = e$est,
                    se = if ("sd" %in% names(e)) e$sd else g("se"),
                    lower = if ("q025" %in% names(e)) e$q025 else g("ci.lower"),
                    upper = if ("q975" %in% names(e)) e$q975 else g("ci.upper"), stringsAsFactors = FALSE)
  if (bayes) {
    out$rhat <- g("rhat"); out$ess <- g("ess"); out$p_incl <- g("p_incl"); out$BF10 <- g("BF10")
    if (any(!is.na(g("q025_adj")))) { out$lower_adj <- g("q025_adj"); out$upper_adj <- g("q975_adj") }
    out$prior <- as.character(g("prior")); out$prior[is.na(out$prior)] <- ""
  } else { out$z <- g("z"); out$pvalue <- g("pvalue") }
  fx <- if ("fixed" %in% names(e)) e$fixed else is.na(out$se)
  out$fixed <- !is.na(fx) & fx
  for (n in c("rhat", "ess", "p_incl", "BF10", "lower_adj", "upper_adj")) if (n %in% names(out) && all(is.na(out[[n]]))) out[[n]] <- NULL
  out
}

partable <- function(e, bayes, title, sections = NULL, spec = NULL) {
  p <- btable(std_partable(e, bayes), title, "birt_partable")
  attr(p, "bayes") <- bayes; attr(p, "sections") <- sections; attr(p, "notes") <- section_notes(p, spec)
  p
}

# notes after a section name, only when they apply: step difficulties for partial credit items,
# the SHASH shape parameters when a distribution has them
section_notes <- function(p, spec) {
  pc <- !is.null(spec$ord$type) && any(unlist(spec$ord$type) != "grm")
  c(Thresholds = if (pc) paste0("partial credit items: step difficulties", if (any(p$op == "|a")) "; |a: step discriminations" else "") else "",
    Distributions = if (any(grepl("^(skew|tail)", p$op))) "SHASH: skew, tail < 1 heavier than normal" else "")
}

# section, group header and row label of each parameter, as lavaan prints them
par_layout <- function(p) {
  self <- p$lhs == p$rhs
  sec <- ifelse(p$op == "=~", "Latent Variables",
         ifelse(p$op == "~" & grepl("^[EV]\\(", p$lhs), "Moderation (E: expectation, V: variance)",
         ifelse(p$op == "~", "Regressions",
         ifelse(p$op == "~~" & !self, "Covariances (correlations)",
         ifelse(p$op == "~~" & self, "Variances",
         ifelse(p$op == "~1", "Intercepts",
         ifelse(p$op %in% c("|", "|a"), "Thresholds",
         ifelse(p$op %in% c("sd", "sd~", "skew", "skew~", "tail"), "Distributions",
                "Other"))))))))
  grp <- ifelse(p$op %in% c("=~", "~", "~~") & !(p$op == "~~" & self), paste(p$lhs, p$op),
         ifelse(sec == "Distributions", p$lhs, ""))
  lab <- ifelse(p$op %in% c("=~", "~"), p$rhs,
         ifelse(p$op == "~~" & !self, p$rhs,
         ifelse(p$op %in% c("|", "|a"), paste0(p$lhs, "|", p$rhs),
         ifelse(grepl("^(sd|skew|tail)", p$op), paste0(p$op, p$rhs), p$lhs))))
  data.frame(section = sec, group = grp, label = lab, stringsAsFactors = FALSE)
}

#' @export
print.birt_partable <- function(x, digits = 3, title = TRUE, ...) {
  if (title && !is.null(attr(x, "title"))) cat(attr(x, "title"), "\n")
  bayes <- isTRUE(attr(x, "bayes")); p <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  L <- par_layout(p)
  num <- if (bayes) intersect(c("est", "se", "lower", "upper", "lower_adj", "upper_adj", "rhat", "ess", "p_incl", "BF10"), names(p))
         else c("est", "se", "z", "pvalue", "lower", "upper")
  lab <- c(est = "Estimate", se = if (bayes) "Post.SD" else "Std.Err", lower = if (bayes) "pi.lower" else "ci.lower",
           upper = if (bayes) "pi.upper" else "ci.upper", lower_adj = "adj.lower", upper_adj = "adj.upper",
           rhat = "Rhat", ess = "ESS", p_incl = "p_incl", BF10 = "BF10", z = "z-value", pvalue = "P(>|z|)")
  fmt <- function(v, n) {
    out <- if (n == "ess") ifelse(is.na(v), "", formatC(round(v), format = "d"))
           else if (n == "BF10") ifelse(is.na(v), "", ifelse(is.infinite(v), "Inf", formatC(v, digits = 2, format = "f")))
           else ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
    out
  }
  cells <- vapply(num, function(n) fmt(p[[n]], n), character(nrow(p)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  cells[p$fixed, setdiff(num, "est")] <- ""
  pri <- if (bayes && "prior" %in% names(p)) ifelse(p$fixed, "", p$prior) else NULL
  rowlab <- ifelse(L$group == "", paste0("    ", L$label), paste0("    ", L$label))
  w0 <- max(nchar(c(rowlab, paste0("  ", L$group))), 20)
  w <- pmax(nchar(lab[num]), apply(cells, 2, function(v) max(nchar(v), 0)))
  line <- function(first, vals, prior = NULL) {
    s <- paste0(sprintf("%-*s", w0, first), paste(sprintf("%*s", w, vals), collapse = "  "))
    if (!is.null(prior) && nzchar(prior)) s <- paste0(s, "  ", prior)
    sub("\\s+$", "", s)
  }
  secs <- unique(L$section); ord <- attr(x, "sections") %||% secs
  for (s in c(intersect(ord, secs), setdiff(secs, ord))) {
    i <- which(L$section == s)
    i <- i[order(match(L$group[i], unique(L$group[i])), grepl(":", L$label[i]))]   # rows of a group together (cross-loadings); main effects first
    nt <- attr(x, "notes")[s]
    cat(sprintf("\n%s%s:\n", s, if (!is.na(nt) && nzchar(nt)) sprintf(" (%s)", nt) else ""))
    cat(line("", lab[num], if (!is.null(pri)) "Prior"), "\n", sep = "")
    last <- NULL
    for (k in i) {
      if (L$group[k] != "" && !identical(L$group[k], last)) { cat("  ", L$group[k], "\n", sep = ""); last <- L$group[k] }
      cat(line(rowlab[k], cells[k, ], if (!is.null(pri)) pri[k]), "\n", sep = "")
    }
  }
  invisible(x)
}

# wide item tables (one row per item; est and se), from the parameter table
item_wide <- function(p, items, factors, binary, ordinal = character()) {
  pick <- function(op, lhs = NULL, rhs = NULL, col) {
    k <- p$op == op & (if (is.null(lhs)) TRUE else p$lhs == lhs) & (if (is.null(rhs)) TRUE else p$rhs == rhs)
    if (any(k)) p[[col]][k][1] else NA_real_
  }
  build <- function(col) {
    out <- data.frame(item = items, type = ifelse(items %in% binary, "binary", ifelse(items %in% ordinal, "ordinal", "continuous")), stringsAsFactors = FALSE)
    for (f in factors) out[[f]] <- vapply(items, function(v) pick("=~", f, v, col), 0)
    out$intercept <- vapply(items, function(v) pick("~1", v, NULL, col), 0)
    out$resid_var <- vapply(items, function(v) pick("~~", v, v, col), 0)
    if (any(p$op == "irt_b")) {
      out$a <- vapply(items, function(v) pick("irt_a", v, NULL, col), 0); out$b <- vapply(items, function(v) pick("irt_b", v, NULL, col), 0)
    }
    nt <- max(c(0, as.integer(sub("^t", "", p$rhs[p$op == "|"]))))
    for (k in seq_len(nt)) out[[sprintf("t%d", k)]] <- vapply(items, function(v) pick("|", v, sprintf("t%d", k), col), 0)
    na <- max(c(0, as.integer(sub("^a", "", p$rhs[p$op == "|a"]))))
    for (k in seq_len(na)[-1]) out[[sprintf("a%d", k)]] <- vapply(items, function(v) pick("|a", v, sprintf("a%d", k), col), 0)
    if (all(is.na(out$intercept))) out$intercept <- NULL
    if (all(is.na(out$resid_var))) out$resid_var <- NULL
    out
  }
  list(est = build("est"), se = build("se"))
}
