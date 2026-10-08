test_that("a small model fits and the marginal likelihood matches brute force", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  d <- sim_rtirt(N = 80, K = 4)
  m <- "theta =~ y1 + y2 + y3 + y4\nspeed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4\nspeed ~ x1"
  fit <- birt(m, d, family = list(speed = ald(0.3)), n_chains = 2, n_iter = 400, n_burn = 200, progress = FALSE, engine = "jags")
  expect_s3_class(fit, "birt")
  expect_equal(nrow(scores(fit)), 80)
  expect_true(all(c("theta", "speed") %in% names(reliability(fit))))
  expect_equal(dim(pointwise_loglik(fit))[2], 80)
  # a draw with moderate loadings: with a logit loading above about 10 the integrand in theta is a
  # step and the Gauss-Hermite rule of the package is off by a few hundredths
  M <- as.matrix(fit$draws); s1 <- which.min(apply(M[, grep("^lam_theta", colnames(M)), drop = FALSE], 1, max))
  P <- birt:::prep_pars(fit); ob <- birt:::obs_matrices(fit)
  LL <- birt:::loglik_compute(P, ob)
  y <- ob$Y[1, ]; t <- ob$T[1, ]; m2 <- birt:::mean_cov(P$F2, s1, 80)[1]
  f <- function(th, s) exp(sum(dbinom(y, 1, plogis(P$d[s1, ] + P$a1[s1, ] * th), log = TRUE)) +
    sum(dnorm(t, P$xi[s1, ] + P$l2[s1, ] * s, P$sig[s1, ], log = TRUE)) + dnorm(th, log = TRUE) +
    birt:::dald_log(s, m2, P$F2$sc[s1], 0.3) - LL[s1, 1])
  inner <- function(th) integrate(Vectorize(function(s) f(th, s)), -4, 4, subdivisions = 2000)$value
  bf <- LL[s1, 1] + log(integrate(Vectorize(inner), -8, 8, subdivisions = 2000)$value)
  expect_equal(LL[s1, 1], bf, tolerance = 1e-3)
  expect_s3_class(ppc(fit, 20), "birt_ppc")
})

test_that("quantile fits print without fit indices", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  d <- sim_rtirt(N = 150, K = 4)
  m <- "theta =~ y1 + y2 + y3 + y4\nspeed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4\nspeed ~ theta + x1"
  q <- birt(m, d, family = list(speed = ald(c(0.25, 0.75))), n_chains = 2, n_iter = 400, n_burn = 200,
            fit_indices = FALSE, progress = FALSE, engine = "jags")
  expect_output(print(q), "Reliability")
  expect_equal(dim(coef(q))[2], 2)
})

test_that("summary shows the priors and the spike-and-slab selection", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  d <- sim_rtirt(N = 200, K = 4, rho = c(0, 0.4, 0, 0.4))
  m <- rtirt_syntax(paste0("y", 1:4), paste0("t", 1:4), cross = "ssp")
  fit <- birt(m, d, n_chains = 2, n_iter = 1000, n_burn = 500, fit_indices = FALSE, progress = FALSE, ssp_p = 0.5, engine = "jags")
  out <- utils::capture.output(print(summary(fit)))
  expect_true(any(grepl("hier(ability)", out, fixed = TRUE)))
  expect_true(any(grepl("$selection", out, fixed = TRUE)))
  expect_true(any(grepl("prior inclusion 0.5", out)))
  e <- estimates(fit); s <- e[!is.na(e$p_incl), ]
  expect_equal(s$BF10, s$p_incl / (1 - s$p_incl))                 # prior odds 1
})

test_that("one chain works (no R-hat), and subset_chains() to one chain", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  d <- sim_rtirt(N = 150, K = 4)
  m <- rtirt_syntax(paste0("y", 1:4), paste0("t", 1:4))
  f1 <- birt(m, d, n_chains = 1, n_iter = 300, n_burn = 150, fit_indices = FALSE, progress = FALSE, engine = "jags")
  expect_true(all(is.na(f1$summary_nodes$rhat)))
  expect_output(print(summary(f1)), "[$]parameters")
  f2 <- birt(m, d, n_chains = 2, n_iter = 300, n_burn = 150, fit_indices = FALSE, progress = FALSE, engine = "jags")
  expect_s3_class(subset_chains(f2, 1), "birt")
})

test_that("a JAGS fit leaves the user's workspace alone (R2jags writes its data to envir)", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  assign("N", "user object", envir = globalenv()); on.exit(rm("N", envir = globalenv()))
  d <- sim_rtirt(N = 100, K = 3)
  birt(rtirt_syntax(paste0("y", 1:3), paste0("t", 1:3)), d, n_chains = 2, n_iter = 200, n_burn = 100, fit_indices = FALSE, progress = FALSE, engine = "jags")
  expect_identical(get("N", envir = globalenv()), "user object")
})

test_that("check_sbc(): the R simulator follows the model, and ranks are computed", {
  d <- sim_hgrm(N = 200, J = 4); d$z <- as.numeric(scale(d$logT)); d$x <- rbinom(200, 1, 0.5)
  m <- "f =~ y1 + y2 + y3 + y4 + x\nE(y1) ~ z\nE(f) ~ z\nV(f) ~ z"
  sp <- birt:::build_spec(m, d, ordered = paste0("y", 1:4), itemtype = c(y1 = "tppcm", y2 = "gpcm")); sp$dp <- dpriors()
  pt <- birt:::build_code(sp)$partab
  set.seed(1); sim <- birt:::sbc_simulate(sp, pt, d)
  expect_true(all(stats::na.omit(unlist(pt$node)) %in% names(sim$truth)))      # every monitored parameter has a true value
  expect_true(all(sim$data$x %in% c(0, 1, NA))); expect_true(all(sim$data$y1 %in% c(sp$ord$levels$y1, NA)))
  expect_equal(is.na(sim$data$y3), is.na(d$y3))
  x <- replicate(4000, birt:::rprior("dt(0, 1, 3) T(0,)")); expect_true(all(x > 0))
  expect_equal(mean(x <= 1), 2 * (pt(1, 3) - 0.5), tolerance = 0.03)          # half-t(3): P(x <= 1)
  expect_equal(sd(replicate(4000, birt:::rprior("dnorm(0, 1/4)"))), 2, tolerance = 0.05)   # precision 1/4: SD 2
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  s <- check_sbc(m, d, n_rep = 2, n_draws = 20, ordered = paste0("y", 1:4), itemtype = c(y1 = "tppcm", y2 = "gpcm"),
                 n_iter = 300, n_burn = 150, progress = FALSE)
  expect_s3_class(s, "birt_sbc"); expect_equal(dim(s$ranks)[1], 2); expect_true(all(s$ranks >= 0 & s$ranks <= 20))
  expect_output(print(s), "Simulation-based calibration")
})

test_that("JAGS fit indices: marginal log-likelihood is exact for narrow posteriors (40 RT items)", {
  skip_on_cran(); skip_if(Sys.which("jags") == "", "JAGS not installed")
  set.seed(7); N <- 40; K <- 40; th <- rnorm(N)
  Tm <- sapply(1:K, function(k) 2 * th + rnorm(N, 0, 0.3)); colnames(Tm) <- paste0("t", 1:K)
  f <- suppressWarnings(birt(paste("f =~", paste(colnames(Tm), collapse = " + ")), as.data.frame(Tm),
                             n_chains = 1, n_iter = 300, n_burn = 150, fit_indices = FALSE, progress = FALSE, engine = "jags"))
  f <- suppressWarnings(add_fit_indices(f, n_draws = 3, progress = FALSE))     # the path birt() uses
  P <- birt:::prep_pars(f, birt:::n_draws_ll(f$spec, 3))
  for (s in seq_len(P$S)) {
    l <- P$l1[s, ]; xi <- P$xi[s, ]; Sig <- outer(l, l) + diag(P$sig[s, ]^2)   # exact normal marginal of 40 RTs
    ex <- vapply(seq_len(N), function(j) { r <- Tm[j, ] - xi
      -0.5 * (K * log(2 * pi) + c(determinant(Sig)$modulus) + sum(r * solve(Sig, r))) }, 0)
    expect_equal(unname(f$loglik[s, ]), ex, tolerance = 1e-6)
  }
})
