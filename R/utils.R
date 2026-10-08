# Internal helpers ---------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# JAGS-safe node suffix for a variable name
safe <- function(x) gsub("[^A-Za-z0-9]", "_", x)

# number as written into JAGS code
num <- function(x) trimws(formatC(x, digits = 10, format = "g"))

# row-wise log-sum-exp of a matrix
lse <- function(M) {
  m <- apply(M, 1, max)
  m + log(rowSums(exp(M - m)))
}

# Asymmetric Laplace distribution AL(mu, sigma, p):
#   f(y) = p (1 - p) / sigma * exp(-rho_p((y - mu) / sigma)),  rho_p(u) = u (p - I(u < 0)).
# Mixture representation (Kozumi & Kobayashi, 2011):
#   y = mu + k1 e + sqrt(k2 sigma e) z,  e ~ Exp(mean sigma), z ~ N(0, 1).
ald_k <- function(p) c(k1 = (1 - 2 * p) / (p * (1 - p)), k2 = 2 / (p * (1 - p)))
ald_var <- function(p) (1 - 2 * p + 2 * p^2) / (p^2 * (1 - p)^2)   # variance / sigma^2
dald_log <- function(y, mu, sigma, p) {
  u <- (y - mu) / sigma
  log(p * (1 - p) / sigma) - u * (p - (u < 0))
}
rald <- function(n, mu, sigma, p) {
  k <- ald_k(p); e <- stats::rexp(n, 1 / sigma)
  mu + k[1] * e + sqrt(k[2] * sigma * e) * stats::rnorm(n)
}

# Empirical reliability of EAP scores: var(EAP) / (var(EAP) + mean PSD^2), which by the
# law of total variance equals var(EAP) / var(latent).
emp_rel <- function(draws) var_rel(colMeans(draws), apply(draws, 2, stats::sd))

# empirical reliability from posterior means m and posterior SDs s: var(m) / (var(m) + mean(s^2))
var_rel <- function(m, s) stats::var(m) / (stats::var(m) + mean(s^2))

# Boxed table printer
# plain table (as print.data.frame, without row names): numbers right-aligned with a fixed
# number of decimals, text columns left-aligned, two spaces between columns
print_table <- function(df, digits = 3) {
  fmt <- function(x) {
    if (!is.numeric(x)) return(ifelse(is.na(x), "", as.character(x)))
    if (is.integer(x)) return(ifelse(is.na(x), "", formatC(x, format = "d", big.mark = "")))
    ifelse(is.na(x), "", ifelse(is.infinite(x), ifelse(x > 0, "Inf", "-Inf"),
                                formatC(x, digits = digits, format = "f")))
  }
  cells <- vapply(df, fmt, character(nrow(df)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  head <- names(df)
  w <- pmax(nchar(head), apply(cells, 2, function(x) max(nchar(x), 0)))
  left <- !vapply(df, is.numeric, TRUE)
  row <- function(x) sub("\\s+$", "", paste(ifelse(left, sprintf("%-*s", w, x), sprintf("%*s", w, x)), collapse = "  "))
  cat(c(row(head), apply(cells, 1, row)), sep = "\n")
}


# max ignoring NA; NA when there is nothing (e.g. R-hat of a single chain)
max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)

pkg_version <- function() tryCatch(as.character(utils::packageVersion("birt")), error = function(e) "")


# Set the RNG seed for the calling function and restore the user's RNG state when it exits
# (so a fixed default seed does not change the global random number stream).
local_seed <- function(seed, envir = parent.frame()) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  restore <- function() {
    if (had) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  }
  do.call(base::on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
  set.seed(seed)
}
