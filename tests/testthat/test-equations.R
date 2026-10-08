test_that("equations() writes ordinal and moderated models in Unicode and LaTeX", {
  skip_if_not_installed("RTMB")
  d <- sim_hgrm(N = 200, J = 5)
  h <- suppressMessages(suppressWarnings(birt("SA =~ y1 + y2 + y3 + y4 + y5\nE(y1 + y2) ~ logT\nE(y1 + y2 + y3 + y4 + y5) ~ pa*logT:SA\nE(SA) ~ logT\nV(SA) ~ logT",
                                              d, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))))
  u <- as.character(equations(h))
  expect_true(any(grepl("logit P(y", u, fixed = TRUE)))
  expect_true(any(grepl("α", u)))                        # alpha (loading moderation)
  expect_true(any(grepl("label pa", u)))
  expect_false(any(grepl("Priors", u)))                  # maximum likelihood
  l <- as.character(equations(h, "latex"))
  expect_equal(l[1], "\\begin{aligned}"); expect_equal(l[length(l)], "\\end{aligned}")
  expect_true(any(grepl("\\alpha_{i", l, fixed = TRUE)))
  expect_false(any(grepl("[^\\x01-\\x7f]", l, perl = TRUE)))         # LaTeX output is ASCII
})

test_that("equations() of a lavaan-syntax fit", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 400, K = 5)
  f <- suppressWarnings(birt(rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5)), d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  u <- as.character(equations(f))
  expect_true(any(grepl("ξᵢ − speed", u)))             # xi_i - speed_j
  expect_true(any(grepl("Corr(ability", u, fixed = TRUE)))
  expect_output(print(equations(f, "latex")), "mathcal\\{N\\}")
})
