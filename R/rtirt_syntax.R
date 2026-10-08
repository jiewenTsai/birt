#' Model syntax for joint response accuracy and response time models
#'
#' Writes the lavaan-style syntax of the usual RT-IRT models from the item and time
#' column names, to pass to [birt()] or to edit.
#'
#' @param items Names of the binary item columns.
#' @param times Names of the (log) response time columns, in the order of `items`.
#' @param cross Cross-loadings of ability on the response times: `"none"`, `"free"`, or
#'   `"ssp"` (spike-and-slab priors; engine `"jags"`).
#' @param structure Relation of ability and speed: `"cor"` (correlated, van der Linden,
#'   2007), `"path"` (`speed ~ ability`), or `"none"` (independent). The default is `"cor"`
#'   without cross-loadings and `"none"` with them (a common cross-loading is the same as the
#'   correlation or the path, so they are not identified together).
#' @param cov Covariates of the latent regression(s).
#' @param cov_on Which latent variables are regressed on `cov`: `"speed"`, `"ability"` or `"both"`.
#' @param speed_var `"free"` (default; the speed loadings are fixed at -1, so the speed
#'   variance is estimated) or `"fixed"` (`speed ~~ 1*speed`, as in the dissertation's
#'   Chapter 5; usually too restrictive for log seconds).
#' @param ability,speed Names of the latent variables.
#' @return The syntax (a character string of class `birt_syntax`, printed with `cat`).
#' @examples
#' items <- paste0("y", 1:5); times <- paste0("t", 1:5)
#' rtirt_syntax(items, times)                                    # van der Linden (2007)
#' rtirt_syntax(items, times, cross = "ssp")                     # cross-loadings, ssp priors
#' rtirt_syntax(items, times, structure = "path", cov = c("x1", "x2"))
#' @export
rtirt_syntax <- function(items, times, cross = c("none", "free", "ssp"), structure = NULL, cov = NULL,
                         cov_on = c("speed", "ability", "both"), speed_var = c("free", "fixed"),
                         ability = "ability", speed = "speed") {
  cross <- match.arg(cross); cov_on <- match.arg(cov_on); speed_var <- match.arg(speed_var)
  if (!is.character(items) || !is.character(times) || !length(items)) stop("items and times must be column names")
  if (length(items) != length(times)) stop(sprintf("%d items but %d times", length(items), length(times)))
  if (is.null(structure)) structure <- if (cross == "none") "cor" else "none"
  structure <- match.arg(structure, c("cor", "path", "none"))
  if (cross != "none" && structure != "none")
    stop("cross-loadings and '", structure, "' are not identified together (a common cross-loading equals the ",
         if (structure == "cor") "correlation" else "path", "); use structure = \"none\"", call. = FALSE)
  wrap <- function(lhs, op, terms) {
    lines <- split(terms, ceiling(seq_along(terms) / 6))
    paste0(lhs, " ", op, " ", paste(vapply(lines, paste, "", collapse = " + "), collapse = " +\n         "))
  }
  out <- c(wrap(ability, "=~", items), wrap(speed, "=~", paste0("-1*", times)))
  if (cross != "none") out <- c(out, wrap(ability, "=~", if (cross == "ssp") sprintf('prior("ssp")*%s', times) else times))
  if (speed_var == "fixed") out <- c(out, sprintf("%s ~~ 1*%s", speed, speed))
  if (structure == "cor") out <- c(out, sprintf("%s ~~ %s", ability, speed))
  if (length(cov) || structure == "path") {
    sp_terms <- c(if (structure == "path") ability, if (cov_on %in% c("speed", "both")) cov)
    if (length(sp_terms)) out <- c(out, wrap(speed, "~", sp_terms))
    if (length(cov) && cov_on %in% c("ability", "both")) out <- c(out, wrap(ability, "~", cov))
  }
  structure(paste0(paste(out, collapse = "\n"), "\n"), class = c("birt_syntax", "character"))
}

#' @export
print.birt_syntax <- function(x, ...) { cat(x); invisible(x) }
