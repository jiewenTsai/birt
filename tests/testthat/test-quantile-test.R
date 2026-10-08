test_that("quantile_test(): Wald test of the scale moderation by the predictor; sup-t bands", {
  skip_if_not_installed("RTMB")
  m <- paste("ability =~ y1 + y2 + y3 + y4 + y5", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5", "ability =~ t1", "V(t1) ~ ability", sep = "\n")
  d <- sim_rtirt(N = 500, K = 5, seed = 3)
  f <- suppressWarnings(birt(m, d, engine = "rtmb", family = list(t1 = shash()), progress = FALSE, control = rtmb_control(method = "aghq")))
  e <- f$estimates; k <- which(e$lhs == "V(t1)" & e$rhs == "ability")
  qt <- quantile_test(f, c(.1, .5, .9)); tt <- qt$test
  r <- tt$target == "t1" & tt$predictor == "ability"
  expect_equal(tt$chisq[r], (e$est[k] / e$se[k])^2, tolerance = 1e-8)        # the test of the curve = the test of V(t1) ~ ability
  expect_equal(tt$df[r], 1L)
  expect_true(all(tt$note[tt$predictor == "speed"] == "constant in p by the model"))
  q <- quantiles(f, c(.1, .25, .5, .75, .9)); a <- q[q$target == "t1" & q$predictor == "ability", ]
  cv <- (a$band.upper - a$est) / a$se
  expect_true(all(cv > 1.96 & cv < stats::qnorm(1 - 0.025 / 5)))           # between pointwise and Bonferroni
  ok <- !is.na(q$se); expect_equal(sqrt(diag(attr(q, "vcov")))[ok], q$se[ok], tolerance = 1e-10, ignore_attr = TRUE)
  expect_error(quantile_test(f, 0.5), "two")
})

test_that("quantile_score(): influence functions of quantile effects", {
  skip_if_not_installed("RTMB"); skip_if_not_installed("strucchange")
  m <- paste("ability =~ y1 + y2 + y3 + y4 + y5", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5", "ability =~ t1", "V(t1) ~ ability", sep = "\n")
  d <- sim_rtirt(N = 500, K = 5, seed = 3)
  fit <- function(dd) suppressWarnings(birt(m, dd, engine = "rtmb", family = list(t1 = shash()), progress = FALSE, control = rtmb_control(method = "aghq")))
  f <- fit(d); P <- c(.1, .9)
  qs <- quantile_score(f, P); e <- attr(qs, "effects"); q <- quantiles(f, P)
  expect_s3_class(qs, "birt_qscore"); expect_equal(nrow(qs), 500)
  expect_true(all(abs(colMeans(qs)) < 1e-4))                                   # sums to zero at the estimates
  expect_equal(e$se, unname(q[colnames(qs), "se"]), tolerance = 1e-10)                  # model-based SE = delta method of quantiles()
  expect_true(all(abs(e$se.robust / e$se - 1) < 0.15))                          # correct model: sandwich close to model-based
  k <- "t1|ability|0.9"; j <- which.max(abs(qs[, k]))
  qj <- quantiles(fit(d[-j, ]), P)
  expect_equal(unname(qj[k, "est"] - q[k, "est"]), unname(-qs[j, k] / 500), tolerance = 0.15)   # case deletion = -IF / N (first order)
  st <- score_test(qs, factor(rep(1:2, 250)))
  expect_s3_class(st, "birt_score_test"); expect_equal(nrow(st), 2)
  expect_error(score_test(quantile_score(f, P, "unconditional"), d$y1), "conditional")
})
