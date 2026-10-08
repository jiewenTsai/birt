test_that("rt_information(): the bands partition the RT information; censoring matches Monte Carlo", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 600, K = 6, seed = 21)
  m <- paste("ability =~ y1 + y2 + y3 + y4 + y5 + y6", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5 + -1*t6",
             "ability =~ t1 + t2", "V(t1) ~ ability", sep = "\n")
  f <- suppressWarnings(birt(m, d, engine = "rtmb", family = list(t1 = shash(), t2 = shash()), progress = FALSE, control = rtmb_control(method = "aghq")))
  ri <- rt_information(f, theta = c(-1, 1), limits = list(slow = .9, fast = .05))
  cr <- cond_reliability(f, theta = c(-1, 1), se = FALSE)
  expect_equal(ri$check$info_partition + ri$check$speed_prior, cr$info_rt, tolerance = 1e-5)
  expect_equal(as.vector(tapply(ri$bands$share, ri$bands$theta, sum)), c(1, 1), tolerance = 1e-12)
  expect_true(all(ri$limits$kept <= 1 + 1e-8 & ri$limits$kept > 0.5))
  expect_true(all(ri$limits$rel <= ri$limits$rel_full + 1e-8))
  # Monte Carlo of the censored information (time limit at u = .9) at theta = 1
  S0 <- f$S; p <- f$rtmb$par; st <- birt:::cr_setup(p, S0, 1, colMeans(S0$X)); set.seed(3); R <- 1e5; I <- matrix(0, 2, 2)
  for (j in seq_along(S0$C)) {
    s <- match(j, S0$sh); dl <- if (is.na(s)) 1 else exp(p$ldl[s]); mu0 <- st$P0$mu[, j]; ls0 <- st$P0$lsc[, j]; e0 <- st$P0$e[, j]
    Tz <- function(z, e) if (is.na(s)) z else birt:::T_std(z, e, dl)
    tt <- mu0 + exp(ls0) * Tz(rnorm(R), e0); cl <- mu0 + exp(ls0) * Tz(qnorm(.9), e0); cens <- tt > cl
    ll <- function(mu, ls, e) { x <- (tt - mu) / exp(ls); y <- (cl - mu) / exp(ls)
      if (is.na(s)) return(ifelse(cens, pnorm(y, lower.tail = FALSE, log.p = TRUE), dnorm(x, log = TRUE) - ls))
      a <- birt:::eff_a(e, dl)
      ifelse(cens, pnorm(sinh(dl * (asinh(y * cosh(a) / dl + sinh(a)) - a)), lower.tail = FALSE, log.p = TRUE), birt:::ld_std(x, e, dl) - ls) }
    h <- 1e-5
    sc <- sapply(list(st$D1, st$D2), function(D) (ll(mu0 + h * D$mu[, j], ls0 + h * D$lsc[, j], e0 + h * D$e[, j]) -
                                                  ll(mu0 - h * D$mu[, j], ls0 - h * D$lsc[, j], e0 - h * D$e[, j])) / (2 * h))
    I <- I + crossprod(sc) / R
  }
  mcinfo <- birt:::cr_rel(list(i11 = I[1, 1], i12 = I[1, 2], i22 = I[2, 2]), st$nF, st$s1)$info
  expect_equal(ri$limits$info_rt[ri$limits$theta == 1 & ri$limits$side == "slow"], mcinfo, tolerance = 0.02)
})
