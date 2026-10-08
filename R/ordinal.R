# Ordinal indicators: graded response (cumulative logits) and partial credit (adjacent-category
# logits) models ---------------------------------------------------------------------------------
#   grm:   P(y > k) = logistic{a (f - delta - b_k)},                 b_1 < b_2 < ...
#   gpcm:  log P(y = k + 1) / P(y = k) = a (f - delta - b_k)          (Muraki, 1992)
#   tppcm: log P(y = k + 1) / P(y = k) = a_k (f - delta - b_k),  a_k = a r_k, r_1 = 1   (Yu, 1991)
# with a the loading of the item (its first step discrimination), delta the shift E(y) ~ z.

ord_types <- c("grm", "gpcm", "tppcm")

# category probabilities (N x K) at values th of (f - delta), for discriminations a (N-vector or
# scalar), thresholds / step difficulties b (K - 1) and step ratios r (K - 1; 1 except tppcm)
# header line for the ordinal indicators of a spec, e.g. "5 ordinal: grm (categories 4)"
ord_desc <- function(sp) sprintf("%d ordinal: %s (categories %s)", length(sp$ordinal), paste(sort(unique(sp$ord$type)), collapse = "/"),
                                 paste(sort(unique(sp$ord$K)), collapse = "/"))

cat_probs <- function(type, a, th, b, r = rep(1, length(b))) {
  K <- length(b) + 1
  if (type == "grm") {
    up <- cbind(1, stats::plogis(a * outer(th, b, "-")), 0)
    return(pmax(up[, -(K + 1), drop = FALSE] - up[, -1, drop = FALSE], 0))
  }
  lp <- if (K == 2) cbind(0, a * r[1] * (th - b[1])) else cbind(0, t(apply(sweep(a * outer(th, b, "-"), 2, r, "*"), 1, cumsum)))
  lp <- lp - apply(lp, 1, max)
  e <- exp(lp); e / rowSums(e)
}

# item information at values th of (f - delta): Var of the score function. grm: Samejima's formula;
# partial credit: the variance of the cumulative step discriminations A_y = sum_{h <= y} a r_h
cat_info <- function(type, a, th, b, r = rep(1, length(b))) {
  P <- cat_probs(type, a, th, b, r)
  if (type == "grm") {
    up <- cbind(1, stats::plogis(a * outer(th, b, "-")), 0); w <- up * (1 - up)
    return(a^2 * rowSums((w[, -ncol(w)] - w[, -1])^2 / pmax(P, 1e-12)))
  }
  A <- matrix(a, length(th), length(b) + 1) * matrix(c(0, cumsum(r)), length(th), length(b) + 1, byrow = TRUE)
  rowSums(P * A^2) - rowSums(P * A)^2
}
