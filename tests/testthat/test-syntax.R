d <- sim_rtirt(N = 50, K = 4)
base <- "theta =~ y1 + y2 + y3 + y4\nspeed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4"

test_that("binary indicators are detected and the code is generated", {
  sp <- birt:::build_spec(paste(base, "theta =~ prior(\"ssp\")*t1 + t2", sep = "\n"), d)
  expect_equal(sp$binary, paste0("y", 1:4))
  expect_equal(sp$cont, paste0("t", 1:4))
  cg <- birt:::build_code(sp)
  expect_match(cg$code, "incl_lam_theta_t1 ~ dbern(p_incl_load)", fixed = TRUE)
  expect_match(cg$code, "lam_speed_t1 <- -1", fixed = TRUE)
  expect_true(all(c("theta", "speed", "sd_speed") %in% cg$monitor))
})

test_that("identification checks", {
  expect_error(birt:::build_spec(paste(base, "theta =~ r*t1 + r*t2 + r*t3 + r*t4\nspeed ~ theta", sep = "\n"), d),
               "not identified")
  expect_silent(birt:::build_spec(paste(base, "theta =~ t1 + t2\nspeed ~ theta", sep = "\n"), d))
  expect_error(birt:::build_spec(paste(base, "theta ~~ speed", sep = "\n"), d, family = list(speed = ald(0.5))),
               "cannot be combined")
  expect_error(birt:::build_spec("a =~ y1\nb =~ y2\nc =~ y3", d), "at most two")
})

test_that("family names resolve by role", {
  sp <- birt:::build_spec(base, d, family = list(speed = ald(0.25), .continuous = ald(0.5), t2 = normal()))
  expect_equal(sp$family$factor$speed$family, "ald")
  expect_equal(sp$family$ind$t1$family, "ald")
  expect_equal(sp$family$ind$t2$family, "normal")
  expect_error(birt:::build_spec(base, d, family = list(y1 = ald(0.5))), "logit link")
  expect_error(birt:::build_spec(base, d, family = list(zz = ald(0.5))), "not in the model")
  expect_error(ald(1.2))
})

test_that("ALD density", {
  f <- function(y) exp(birt:::dald_log(y, 0.3, 0.7, 0.2))
  expect_equal(integrate(f, -Inf, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(f, -Inf, 0.3)$value, 0.2, tolerance = 1e-6)
})

test_that("default priors (dpriors) and a fixed inclusion probability (ssp_p)", {
  sp <- birt:::build_spec(paste(base, "theta =~ prior(\"ssp\")*t1 + t2", sep = "\n"), d)
  sp$dp <- dpriors(cross = "dnorm(0, 100)"); sp$ssp_p <- 0.3
  cg <- birt:::build_code(sp)
  expect_match(cg$code, "lam_theta_t2 ~ dnorm(0, 100)", fixed = TRUE)
  expect_match(cg$code, "p_incl_load <- 0.3", fixed = TRUE)
  expect_false("p_incl_load" %in% cg$monitor)
  pt <- cg$partab
  expect_equal(pt$prior[pt$lhs == "theta" & pt$rhs == "y1"], "hier(theta)")
  expect_match(cg$code, "lam_theta_y1 ~ dlnorm(mu_load_theta, 1)", fixed = TRUE)
  expect_true("mu_load_theta" %in% cg$monitor)
  sp$dp <- dpriors(loading = "dlnorm(0, 1)")
  expect_match(birt:::build_code(sp)$code, "lam_theta_y1 ~ dlnorm(0, 1)", fixed = TRUE)
  expect_error(dpriors(cross = "hier"), "JAGS distributions")
  expect_equal(pt$prior[pt$lhs == "theta" & pt$rhs == "t1"], "ssp")
  expect_equal(pt$prior[pt$lhs == "y1" & pt$op == "~1"], "dnorm(0, 1/9)")
  expect_error(dpriors(beta = "normal(0, 1)"), "JAGS distributions")
})

test_that("defined parameters and constraints are refused, not ignored", {
  expect_error(birt:::parse_model("f =~ y1 + a*y2 + y3\nb := a^2"), "not supported")
  expect_error(birt:::parse_model("f =~ y1 + a*y2 + y3\na > 0"), "not supported")
})

test_that("a label times a number is refused, not silently read as plain equality", {
  for (m in c("f =~ y1 + (p*0.025)*y2", "f =~ y1 + 0.5*p*y2", "f =~ y1 + p*0.5*y2"))
    expect_error(parse_model(m), "one value or one label")
  expect_silent(parse_model("f =~ y1 + p*y2 + p*y3 + prior(\"dnorm(0, 1)\")*y4 + start(0.5)*y5 + -1*y6 +\n  0*y7"))
})

test_that("shared labels: loadings with loadings, regressions with regressions, nothing else", {
  expect_error(parse_model("f =~ p*y1 + y2\ng =~ t1 + t2\nf ~ p*x"), "loading .* and a regression")
  expect_error(parse_model("f =~ y1 + y2\ng =~ t1 + t2\nf ~~ c*g\nf ~~ c*f"), "variance or covariance")
  sp <- parse_model("f =~ y1 + y2 + y3\nf ~ b*x + b*z")
  d <- data.frame(y1 = c(0, 1, 1, 0), y2 = c(1, 0, 1, 0), y3 = c(1, 1, 0, 0), x = 1:4, z = 4:1)
  code <- build_code(build_spec("f =~ y1 + y2 + y3\nf ~ b*x + b*z", d))
  expect_match(code$code, "beta_f_x <- labb_b", fixed = TRUE)
  expect_match(code$code, "beta_f_z <- labb_b", fixed = TRUE)
  expect_equal(sum(grepl("labb_b ~", strsplit(code$code, "\n")[[1]], fixed = TRUE)), 1)
  expect_equal(code$partab$node[code$partab$op == "~"], c("labb_b", "labb_b"))
})

test_that("starting values follow the data and positive loadings start positive", {
  set.seed(3); N <- 400; f <- rnorm(N)
  d <- data.frame(sapply(1:6, function(i) 4 + 0.2 * f + rnorm(N, 0, 0.4)))
  names(d) <- paste0("t", 1:6)
  sp <- birt:::build_spec("speed =~ t1 + t2 + t3 + t4 + t5 + t6", d)
  sv <- birt:::start_values(sp, d)
  expect_true(all(abs(sv$load$speed / 0.2 - 1) < 0.3))
  expect_gt(cor(sv$person$speed, f), 0.7)          # alpha about .6 with six items
  cg <- birt:::build_code(sp); mk <- birt:::make_inits(sp, d, cg, 1)
  expect_true(all(vapply(1:50, function(ch) { i <- mk(ch); min(unlist(i[grep("^lam_", names(i))])) > 0 }, TRUE)))
})

test_that("syntax checks: label + prior, V() labels, 0*z, shared labels, all anchors", {
  m <- "theta =~ y1 + y2 + y3 + y4 + y5"
  d <- sim_hgrm(N = 150, J = 5); d$z <- as.numeric(scale(d$logT))
  expect_error(birt:::parse_model(paste(m, 'theta =~ prior("dnorm(0, 1)")*a1*y1', sep = "\n")), "label and prior")
  expect_error(birt:::build_spec(paste(m, "V(theta) ~ v1*z", sep = "\n"), d, ordered = TRUE), "V")
  expect_null(birt:::build_spec(paste(m, "V(theta) ~ 0*z", sep = "\n"), d, ordered = TRUE)$family$factor$theta$scale)
  expect_error(birt:::parse_model(paste(m, "theta ~ c1*x", "E(y1) ~ c1*x", sep = "\n")), "shared")
  expect_error(birt:::hgrm_syntax(m, "z", "free", "free", "none", paste0("y", 1:5)), "anchor")
})

test_that("ordered = TRUE takes integer indicators with 3-10 categories; itemtype needs ordered", {
  d <- sim_hgrm(N = 150, J = 3); d$cnt <- rpois(150, 30)
  expect_message(s <- birt:::build_spec("theta =~ y1 + y2 + y3 + cnt", d, ordered = TRUE), "10")
  expect_equal(s$cont, "cnt"); expect_equal(s$ordinal, c("y1", "y2", "y3"))
  expect_error(birt:::build_spec("theta =~ y1 + y2 + y3", d, itemtype = "gpcm"), "ordered")
})

test_that("JAGS truncates a user prior on a measurement loading at 0", {
  d <- sim_hgrm(N = 150, J = 3)
  s <- birt:::build_spec('theta =~ prior("dnorm(0, 1)")*y1 + y2 + y3', d); s$dp <- dpriors()
  expect_match(birt:::build_code(s)$code, "T(0,)", fixed = TRUE)
})
