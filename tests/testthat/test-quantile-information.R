test_that("quantile_information(): normal location effect gives f(q)^2 beta^2 / (p (1 - p)) with beta the cross-loading", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 400, K = 4, seed = 1)
  m <- paste("theta =~ y1 + y2 + y3 + y4", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4", "theta =~ t1 + t2 + t3 + t4", sep = "\n")
  f <- birt(m, d, control = rtmb_control(method = "aghq"), progress = FALSE)
  qi <- quantile_information(f, p = c(.1, .5, .9), theta = c(-1, 1))
  est <- f$estimates; lam <- est$est[est$lhs == "theta" & est$op == "=~" & est$rhs %in% f$S$C]
  it <- qi$items
  expect_equal(it$beta, lam[match(it$item, f$S$C)], tolerance = 1e-6)
  sd <- sqrt(est$est[est$op == "~~" & est$lhs == est$rhs & est$lhs %in% f$S$C])[match(it$item, f$S$C)]   # residual variances
  expect_equal(it$density, dnorm(qnorm(it$p)) / sd, tolerance = 1e-8)
  expect_equal(it$info, dnorm(qnorm(it$p))^2 / (it$p * (1 - it$p)) * lam[match(it$item, f$S$C)]^2 / sd^2, tolerance = 1e-6)
  expect_true(all(qi$curve$efficiency > 0 & qi$curve$efficiency < 1))
  expect_equal(qi$best$p, c(.5, .5))
  expect_s3_class(qi, "birt_qinfo")
  expect_output(print(qi), "Best p")
})

test_that("quantile_information(): SHASH with a scale effect; beta and f(q) match the quantile function", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 600, K = 6, seed = 21)
  m <- paste("ability =~ y1 + y2 + y3 + y4 + y5 + y6", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5 + -1*t6",
             "ability =~ t1 + t2", "V(t1) ~ ability", sep = "\n")
  f <- suppressWarnings(birt(m, d, family = list(t1 = shash(), t2 = shash()), progress = FALSE, control = rtmb_control(method = "aghq")))
  pp <- c(.1, .5, .9); th <- 1
  qi <- quantile_information(f, p = pp, theta = th)
  S0 <- f$S; p <- f$rtmb$par; xb <- colMeans(S0$X); h <- 1e-4
  st <- lapply(th + c(-h, 0, h), function(t) birt:::cr_setup(p, S0, t, xb))
  for (j in 1:2) {
    dl <- exp(p$ldl[match(j, S0$sh)])
    qf <- function(s, u) s$P0$mu[, j] + exp(s$P0$lsc[, j]) * birt:::T_std(qnorm(u), s$P0$e[, j], dl)
    r <- qi$items[qi$items$item == S0$C[j], ]
    expect_equal(r$beta, (qf(st[[3]], pp) - qf(st[[1]], pp)) / (2 * h), tolerance = 1e-5)
    expect_equal(r$density, 1 / ((qf(st[[2]], pp + 1e-6) - qf(st[[2]], pp - 1e-6)) / 2e-6), tolerance = 1e-5)
    expect_equal(r$info, r$density^2 * r$beta^2 / (pp * (1 - pp)), tolerance = 1e-5)
  }
  # the scale effect of t1 makes its quantiles move unequally
  b1 <- qi$items$beta[qi$items$item == "t1"]; expect_true(b1[1] > b1[3])
  expect_true(all(qi$curve$info <= qi$curve$info_rt))
  expect_error(quantile_information(f, p = 1), "quantile levels")
  expect_equal(sum(quantile_information(f, theta = 0)$curve$p %in% c(0.75, 0.9)), 2)
})

test_that("sim_rtirt(scale = 0) is the default simulation; scale moves the log residual SD with ability", {
  expect_identical(sim_rtirt(N = 50, K = 3, seed = 3), sim_rtirt(N = 50, K = 3, seed = 3, scale = 0))
  expect_false(identical(sim_rtirt(N = 50, K = 3, seed = 3), sim_rtirt(N = 50, K = 3, seed = 3, scale = 0.3)))
})
