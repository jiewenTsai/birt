m_cross <- "
theta =~ y1 + y2 + y3 + y4 + y5
speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5
theta =~ t1 + t2 + t3 + t4 + t5
"

test_that("AGHQ marginal likelihood matches brute-force integration", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 200, K = 5, seed = 3)
  d$y2[2] <- NA; d$t3[5] <- NA                                     # missing indicators
  fam <- list(.continuous = shash(), t2 = normal(), t3 = shash(skew = "theta"))
  f <- birt(paste(m_cross, "V(t2 + t3) ~ theta", sep = "\n"), d, family = fam, engine = "rtmb", control = rtmb_control(k = 11), progress = FALSE)
  expect_equal(f$rtmb$convergence, 0)
  S <- f$S; p <- f$rtmb$par
  g <- seq(-7, 7, length.out = 241); h <- diff(g[1:2]); G <- expand.grid(z1 = g, z2 = g)
  total <- 0
  for (i in seq_len(S$N)) {
    rows <- rep(i, nrow(G))
    ll <- ind_ll(p, S, lat_values(p, S, G$z1, G$z2, rows), rows) + dnorm(G$z1, log = TRUE) + dnorm(G$z2, log = TRUE)
    total <- total + max(ll) + log(sum(exp(ll - max(ll))) * h^2)
  }
  expect_equal(-f$rtmb$objective, total, tolerance = 1e-4)
})

test_that("the normal model recovers the cross-loadings; quantile effects are constant", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 600, K = 5, rho = c(0, 0.3, 0, 0.3, 0.3), seed = 4)
  f <- birt(m_cross, d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
  expect_s3_class(f, "birt_rtmb")
  e <- estimates(f); cl <- e[e$op == "=~" & e$lhs == "theta" & grepl("^t", e$rhs), ]
  expect_true(all(abs(cl$est - -c(0, 0.3, 0, 0.3, 0.3)) < 3 * cl$se))
  q <- quantiles(f, c(0.1, 0.9)); qt <- q[q$predictor == "theta", ]
  expect_equal(qt$est[qt$p == 0.1], qt$est[qt$p == 0.9])
  expect_equal(qt$est[qt$p == 0.1], cl$est); expect_equal(qt$se[qt$p == 0.1], cl$se)
  expect_equal(as.numeric(logLik(f)), -f$fit_indices$Deviance / 2)
  expect_equal(AIC(f), f$fit_indices$AIC)
  expect_equal(nrow(scores(f)), 600); expect_true(all(reliability(f) > 0 & reliability(f) < 1))
  expect_output(print(summary(f)), "Latent Variables:")
  expect_true(is.numeric(coef(f)) && "theta=~y1" %in% names(coef(f)))
  expect_output(convergence(f), "max_abs_gradient")
  expect_message(plot(f), "no quantile effects")
})

test_that("SHASH quantile effects equal numerical derivatives of the conditional quantiles", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 800, K = 5, seed = 5)
  m <- "theta =~ y1 + y2 + y3 + y4 + y5\nspeed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5\nspeed ~ theta + x1\nV(speed) ~ theta + x1"
  # (normal data: the skewness of speed is weakly identified and may sit at its bound)
  f <- suppressWarnings(birt(m, d, family = list(speed = shash(skew = "theta")), engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  S <- f$S; p <- f$rtmb$par; pq <- 0.8
  # q_p(speed | theta, x1) = gam theta + b x1 + s exp(a1 theta + a2 x1) T(z_p; eps + eta theta, delta)
  qf <- function(th, x) p$gam * th + p$beta2[1] * x +
    exp(p$lsd2 + p$alpha2[1] * th + p$alpha2[2] * x) * T_std(qnorm(pq), p$eps2 + p$eta2[1] * th, exp(p$ldl2))
  gq <- statmod::gauss.quad.prob(15, "normal"); h <- 1e-5
  ame_th <- mean(sapply(d$x1, function(x) sum(gq$weights * (qf(gq$nodes + h, x) - qf(gq$nodes - h, x)) / (2 * h))))
  ame_x <- mean(sapply(d$x1, function(x) sum(gq$weights * (qf(gq$nodes, x + h) - qf(gq$nodes, x - h)) / (2 * h))))
  q <- quantiles(f, pq)
  expect_equal(q$est[q$target == "speed" & q$predictor == "theta"], ame_th, tolerance = 1e-6)
  expect_equal(q$est[q$target == "speed" & q$predictor == "x1"], ame_x, tolerance = 1e-6)
  expect_true(all(q$se[q$target == "speed"] > 0))
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(f, targets = "speed"))
})

test_that("RTMB engine refuses what it cannot fit", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 50, K = 4)
  m2 <- "theta =~ y1 + y2 + y3 + y4\nspeed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4\nspeed ~ x1"
  expect_error(birt(m2, d, family = list(speed = ald(0.3)), engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "shash")
  expect_error(birt(m2, d, family = list(speed = ald(c(0.3, 0.5))), engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "one model for all quantiles")
  expect_error(birt(m2, d, family = list(speed = shash()), engine = "jags"), "engine = \"rtmb\"")
  expect_error(birt(m2, d, family = list(theta = shash()), engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "theta is normal")
})

test_that("rtirt_syntax writes the usual models", {
  it <- paste0("y", 1:3); tm <- paste0("t", 1:3)
  s <- rtirt_syntax(it, tm)
  expect_match(s, "ability =~ y1 + y2 + y3", fixed = TRUE); expect_match(s, "speed =~ -1*t1", fixed = TRUE)
  expect_match(s, "ability ~~ speed", fixed = TRUE)
  expect_match(rtirt_syntax(it, tm, cross = "ssp"), 'prior("ssp")*t2', fixed = TRUE)
  expect_match(rtirt_syntax(it, tm, structure = "path", cov = "x1"), "speed ~ ability + x1", fixed = TRUE)
  expect_match(rtirt_syntax(it, tm, cross = "free", speed_var = "fixed"), "speed ~~ 1*speed", fixed = TRUE)
  expect_error(rtirt_syntax(it, tm, cross = "free", structure = "cor"), "not identified")
  expect_error(rtirt_syntax(it, tm[1:2]), "3 items but 2 times")
  expect_type(parse_model(s), "list")
})

test_that("data checks: 1/2 coding, raw seconds, log_rt and id", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 300, K = 5, seed = 6)
  m <- rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5))
  d2 <- d; d2$y1 <- d2$y1 + 1
  expect_error(birt(m, d2, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "not coded 0/1")
  d3 <- d; for (v in paste0("t", 1:5)) d3[[v]] <- exp(d3[[v]])
  expect_warning(build_spec(m, d3), "raw response times")
  f <- birt(m, d3, engine = "rtmb", log_rt = TRUE, id = 1001:1300, progress = FALSE, control = rtmb_control(method = "aghq"))
  g <- birt(m, d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
  expect_equal(f$fit_indices$Deviance, g$fit_indices$Deviance, tolerance = 1e-6)
  expect_equal(scores(f)$id, 1001:1300)
  expect_output(print(f), "Log-transformed: t1, t2, t3, t4, ... (5)", fixed = TRUE)
  d3$t1[1] <- -1
  expect_error(birt(m, d3, engine = "rtmb", log_rt = TRUE, progress = FALSE, control = rtmb_control(method = "aghq")), "positive seconds")
  expect_output(compare(a = f, b = g), "dAIC")
})

test_that("rtmb_code() writes a script that reproduces the fit", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 400, K = 5, seed = 2); d$y2[3] <- NA
  m <- rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5), cross = "free")
  f <- birt(paste(m, "V(t1 + t3 + t4 + t5) ~ ability", "V(t2) ~ x1", sep = "\n"), d, engine = "rtmb", family = list(.continuous = shash(), t2 = normal()), progress = FALSE, control = rtmb_control(method = "aghq"))
  e <- new.env(); e$data <- d
  utils::capture.output(eval(parse(text = rtmb_code(f)), e))
  expect_equal(e$optA$objective, f$rtmb$objective, tolerance = 1e-5)
  S <- f$S; p <- f$rtmb$par; z1 <- rnorm(400); z2 <- rnorm(400)
  expect_equal(sum(e$joint(p, z1, z2, 1:400)),
               sum(ind_ll(p, S, lat_values(p, S, z1, z2, 1:400), 1:400) + dnorm(z1, log = TRUE) + dnorm(z2, log = TRUE)))
  e2 <- new.env(); e2$data <- d
  eval(parse(text = utils::capture.output(rtmb_code(m, data = d, estimate = FALSE))), e2)
  expect_true(is.function(e2$nll))
  # correlated latent variables with a free SD of speed: speed = exp(lsd2) * (r z1 + sqrt(1 - r^2) z2)
  fc <- birt(rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5)), d, progress = FALSE, control = rtmb_control(method = "aghq"))
  ec <- new.env(); ec$data <- d
  utils::capture.output(eval(parse(text = rtmb_code(fc)), ec))
  expect_equal(ec$optA$objective, fc$rtmb$objective, tolerance = 1e-5)
})

test_that("prior strings become normalized log densities", {
  for (s in c("dnorm(0, 1/9)", "dlnorm(0, 4)", "dt(0, 1, 3) T(0,)", "dunif(-1, 1)", "dgamma(2, 3)", "dbeta(2, 3)", "dnorm(1, 4) T(-1, 2)")) {
    p <- prior_density(s)
    expect_equal(stats::integrate(function(x) exp(p$f(x)), max(p$lo, -200), min(p$hi, 200), subdivisions = 1000)$value, 1, tolerance = 1e-4)
  }
  expect_equal(prior_density("dnorm(0, 1/9)")$f(1), stats::dnorm(1, 0, 3, log = TRUE))
  expect_error(prior_density("normal(0, 1)"), "cannot read")
  expect_error(prior_density("dweib(1, 1)"), "not supported")
})

test_that("ELGM: Gaussian outer layer = posterior mode + Laplace; quadrature fit has priors, BF and mixtures", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 250, K = 5, seed = 21)
  m <- rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5))
  g <- birt(m, d, engine = "rtmb", control = rtmb_control(method = "elgm", hyper = character()), progress = FALSE)
  E <- g$rtmb$elgm; r <- g$rtmb
  expect_equal(length(E$prob), 1)
  ok <- !r$ab
  expect_equal(E$lognc, r$logpost + 0.5 * sum(ok) * log(2 * pi) + 0.5 * determinant(r$vcov[ok, ok])$modulus[1], tolerance = 1e-8)
  e <- g$estimates; i <- e$op == "~1" & e$lhs %in% g$S$B
  expect_equal(e$est[i], unname(r$par$d[seq_len(sum(i))]))     # identity parameters: mean = mode
  f <- birt(m, d, engine = "rtmb", control = rtmb_control(method = "elgm"), progress = FALSE)
  expect_equal(f$rtmb$elgm$hyper, c("lsd2", "atr"))
  expect_equal(length(f$rtmb$elgm$prob), 81)
  expect_equal(sum(f$rtmb$elgm$prob), 1)
  ef <- f$estimates
  expect_true(all(ef$prior[!ef$fixed & ef$op %in% c("=~", "~1", "~~")] != ""))
  expect_true(all(ef$ci.lower[!ef$fixed] < ef$est[!ef$fixed] & ef$est[!ef$fixed] < ef$ci.upper[!ef$fixed]))
  rr <- ef[ef$op == "~~" & ef$lhs != ef$rhs, ]
  expect_true(rr$ci.lower > -1 && rr$ci.upper < 1)
  expect_output(print(summary(f)), "Prior")
  expect_output(print(f), "log marginal likelihood")
  expect_output(convergence(f), "outer quadrature")
  expect_true(all(is.finite(scores(f)$speed_psd)))
  m0 <- paste0(m, "ability ~~ 0*speed\n")
  h <- birt(sub("ability ~~ speed\n", "", m, fixed = TRUE), d, engine = "rtmb", control = rtmb_control(method = "elgm"), progress = FALSE)
  expect_output(tab <- compare(cor = f, indep = h), "BF_best")
  expect_equal(tab$model[which.max(tab$logML)], "cor")
  expect_error(compare(f, birt(m, d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))), "marginal likelihoods")
})

test_that("ELGM grid coverage: the outer-mass warning means the same width for every k_hyper", {
  skip_if_not_installed("RTMB")
  # a Gaussian posterior s times wider than the grid assumes, on the k-node rule (1 and 2 directions)
  edge <- function(k, s, nq = 1) {
    gh <- statmod::gauss.quad.prob(k, "normal"); Zq <- as.matrix(expand.grid(rep(list(gh$nodes), nq)))
    lw <- rowSums(log(as.matrix(expand.grid(rep(list(gh$weights), nq))))) + rowSums(dnorm(Zq, 0, s, log = TRUE) - dnorm(Zq, log = TRUE))
    p <- exp(lw - max(lw)); birt:::elgm_edge(p / sum(p), Zq, gh)
  }
  for (k in c(3, 5, 9, 11)) {
    g <- edge(k, 1); expect_equal(g$mass, g$gauss, tolerance = 1e-10); expect_false(g$wide)   # Gaussian: the rule's own weights
    expect_false(edge(k, 1.3)$wide); expect_true(edge(k, 2.5)$wide)
    expect_false(any(edge(k, 1.3, 2)$wide)); expect_true(all(edge(k, 2.5, 2)$wide))
  }
  expect_true(edge(3, 1.6, 5)$wide[1])      # k = 3 on 5 directions (the old total-mass rule could not fire: edge0 = .87)
  # a well-determined posterior: no warning with k_hyper = 3 or 5
  m <- rtirt_syntax(paste0("y", 1:3), paste0("t", 1:3)); d <- sim_rtirt(N = 200, K = 3, seed = 1)
  for (k in c(3, 5)) expect_no_warning(f <- birt(m, d, progress = FALSE, control = rtmb_control(method = "elgm", k_hyper = k)))
  expect_false(any(f$rtmb$elgm$coverage$wide))
  # N = 25: the posterior of the log speed SD has a long tail towards 0 (16% of the mass on the outer nodes with k = 5)
  expect_warning(f <- birt(m, sim_rtirt(N = 25, K = 3, seed = 3), progress = FALSE, control = rtmb_control(method = "elgm", k_hyper = 5)),
                 "too wide for the outer grid along lsd2")
  expect_gt(f$rtmb$elgm$coverage$mass[1], 0.15)
})

test_that("ELGM uses prior() modifiers and one-factor models", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 250, K = 5, rho = c(0, 0.3, 0, 0.3, 0.3), seed = 22)
  m <- 'ability =~ y1 + y2 + y3 + y4 + y5
speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5
ability =~ prior("dnorm(0, 10000)")*t2 + t4'
  f <- suppressMessages(birt(m, d, engine = "rtmb", control = rtmb_control(method = "elgm"), progress = FALSE))
  e <- f$estimates; cl <- e[e$op == "=~" & e$lhs == "ability" & e$rhs %in% c("t2", "t4"), ]
  expect_true(abs(cl$est[cl$rhs == "t2"]) < 0.03 && cl$se[cl$rhs == "t2"] < 0.012)    # prior SD 0.01
  expect_equal(cl$prior, c("dnorm(0,10000)", "dnorm(0, 1)"))
  expect_true(abs(cl$est[cl$rhs == "t4"]) > 0.1)
  expect_message(birt(m, d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "ignores prior")
  expect_error(birt(sub("dnorm(0, 10000)", "ssp", m, fixed = TRUE), d, engine = "rtmb", control = rtmb_control(method = "elgm"), progress = FALSE), "ssp")
  f1 <- birt("ability =~ y1 + y2 + y3 + y4 + y5", d, engine = "rtmb", control = rtmb_control(method = "elgm"), progress = FALSE)
  expect_equal(length(f1$rtmb$elgm$prob), 1)                     # no hyperparameters: Gaussian outer layer
  expect_true(all(f1$reliability > 0))
})

test_that("plot() of an RTMB fit without quantile effects does not fail", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 500, K = 6, seed = 11)
  f <- birt(paste("ability =~", paste0("y", 1:6, collapse = " + ")), d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_message(plot(f), "no quantile effects")
  # SHASH speed without scale or skewness moderation: constant effects, and the message does not say "normal"
  m2 <- "ability =~ y1 + y2 + y3\nspeed =~ -1*t1 + -1*t2 + -1*t3\nspeed ~ ability"
  g <- suppressWarnings(birt(m2, d, family = list(speed = shash()), progress = FALSE, control = rtmb_control(method = "aghq")))
  out <- capture.output(print(summary(g)$parameters))
  expect_true(any(grepl("^Distributions \\(SHASH: skew, tail < 1 heavier than normal\\):", out)))
  expect_message(plot(g), "location shift")
  expect_false(any(grepl("normal residuals", capture.output(plot(g), type = "message"))))
})

test_that("a diverging discrimination is named in a warning", {
  skip_if_not_installed("RTMB")
  set.seed(1); d <- sim_rtirt(N = 200, K = 5)                        # y4 runs off (loading > 800)
  w <- character()
  withCallingHandlers(birt(paste("ability =~", paste0("y", 1:5, collapse = " + ")), d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")),
                      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("diverging discrimination of y4", w)))
})

test_that("the default fit is RTMB with ELGM; maximum likelihood and JAGS are chosen explicitly", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 150, K = 4, seed = 2)
  f <- suppressWarnings(birt("theta =~ y1 + y2 + y3 + y4", d, progress = FALSE))
  expect_s3_class(f, "birt_rtmb"); expect_false(is.null(f$rtmb$elgm))
  expect_error(estfun.birt_rtmb(f), "maximum likelihood")
  expect_identical(rtmb_control()$method, "elgm")
})

test_that("maximum likelihood says that it ignores prior(\"ssp\") on loadings", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 200, K = 4, seed = 2)
  expect_message(birt(rtirt_syntax(paste0("y", 1:4), paste0("t", 1:4), cross = "ssp"), d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")),
                 "ignores prior")
})

test_that("quantile effects include E() moderation of a continuous indicator", {
  skip_if_not_installed("RTMB"); skip_if_not_installed("statmod")
  dr <- sim_rtirt(N = 400, K = 4, seed = 3); set.seed(3); dr$z <- rnorm(400)
  f <- suppressWarnings(birt("theta =~ y1 + y2 + y3 + y4\nspeed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4\nE(t1) ~ z\nE(t1) ~ z:speed",
                             dr, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  q <- quantiles(f, 0.5)
  expect_true(any(q$target == "t1" & q$predictor == "z"))
})

test_that("rtmb_code() / jags_code() from syntax with ordered = and itemtype =; the script reproduces ordinal fits", {
  skip_if_not_installed("RTMB")
  hd <- sim_hgrm(N = 300, J = 4, seed = 2); hd$z <- as.numeric(scale(hd$logT)); hd$y2[5] <- NA
  m <- paste("f =~ y1 + y2 + y3 + y4", "f ~ z", "E(y1 + y2) ~ z", "E(y3) ~ z:f", sep = "\n")
  e0 <- new.env(); e0$data <- hd
  eval(parse(text = utils::capture.output(rtmb_code(m, hd, ordered = TRUE, estimate = FALSE))), e0)
  expect_true(is.function(e0$nll))
  expect_true(any(grepl("dcat|ilogit|logit", utils::capture.output(jags_code(m, hd, ordered = TRUE, itemtype = "gpcm")))))
  for (ty in c("grm", "gpcm", "tppcm")) {
    f <- suppressWarnings(birt(m, hd, ordered = TRUE, itemtype = ty, control = rtmb_control(method = "aghq"), progress = FALSE))
    e <- new.env(); e$data <- hd
    utils::capture.output(eval(parse(text = rtmb_code(f)), e))
    expect_equal(2 * e$optA$objective, 2 * f$rtmb$objective, tolerance = 1e-5)
    S <- f$S; p <- f$rtmb$par; z1 <- rnorm(300)
    expect_equal(sum(e$joint(p, z1, 0, 1:300)), sum(ind_ll(p, S, lat_values(p, S, z1, 0, 1:300), 1:300) + dnorm(z1, log = TRUE)), tolerance = 1e-10)
  }
})

test_that("ELGM hyper = \"auto\" uses at most 3 outer directions and says which are Gaussian", {
  skip_if_not_installed("RTMB")
  set.seed(5); N <- 250; d <- sim_rtirt(N = N, K = 4, seed = 5)
  m <- paste(rtirt_syntax(paste0("y", 1:4), paste0("t", 1:4)), "V(speed) ~ x1", sep = "\n")
  expect_message(f <- suppressWarnings(birt(m, d, family = list(speed = shash(skew = "ability")), progress = FALSE)),
                 "at most 3 hyperparameters \\(eps2, ldl2, eta2\\); alpha2 is treated as Gaussian")
  expect_equal(f$rtmb$elgm$nq, 3); expect_equal(length(f$rtmb$elgm$prob), 125)
})
