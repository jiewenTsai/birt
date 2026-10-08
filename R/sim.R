#' Simulate joint response accuracy and response time data
#'
#' 2PL accuracy and lognormal response times with item-specific cross-loadings of
#' ability on log response time, optionally an effect of ability on the log residual SD of the
#' log response time (scale), and speed depending on two covariates.
#' @param N Persons.
#' @param K Items.
#' @param rho Cross-loadings of ability on log RT (length 1 or K).
#' @param scale Effects of ability on the log residual SD of log RT (length 1 or K; 0: none, as in
#'   `V(t) ~ ability`).
#' @param sd_speed SD of the speed residual.
#' @param beta Effects of covariates `x1`, `x2` on speed.
#' @param seed Random seed.
#' @return A data frame with `y1..yK` (0/1), `t1..tK` (log RT), `x1`, `x2`.
#' @export
sim_rtirt <- function(N = 500, K = 10, rho = rep(c(0, 0.3), length.out = K), scale = 0, sd_speed = 0.3,
                      beta = c(0.1, 0), seed = 1) {
  local_seed(seed)
  rho <- rep_len(rho, K); scale <- rep_len(scale, K)
  x1 <- stats::rnorm(N); x2 <- stats::rnorm(N)
  theta <- stats::rnorm(N); speed <- beta[1] * x1 + beta[2] * x2 + stats::rnorm(N, 0, sd_speed)
  a <- stats::runif(K, 0.8, 2); d <- stats::rnorm(K); xi <- stats::rnorm(K, 4, 0.3); sig <- stats::runif(K, 0.4, 0.6)
  Y <- matrix(stats::rbinom(N * K, 1, stats::plogis(outer(theta, a) + matrix(d, N, K, byrow = TRUE))), N, K)
  Tm <- matrix(xi, N, K, byrow = TRUE) - speed - outer(theta, rho) + matrix(stats::rnorm(N * K), N, K) * exp(log(matrix(sig, N, K, byrow = TRUE)) + outer(theta, scale))
  colnames(Y) <- paste0("y", 1:K); colnames(Tm) <- paste0("t", 1:K)
  data.frame(Y, Tm, x1 = x1, x2 = x2)
}
