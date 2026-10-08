test_that("conditional reliability: 2PL information and the normal RT formula with speed as a nuisance", {
  skip_if_not_installed("RTMB")
  set.seed(1); d <- sim_rtirt(N = 300, K = 4)
  f1 <- birt("theta =~ y1 + y2 + y3 + y4", d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
  th <- c(-2, 0, 1.5); cr <- cond_reliability(f1, theta = th, se = FALSE, bins = 0)
  e <- estimates(f1); a <- e$est[e$op == "=~"]; dd <- e$est[e$op == "~1"]
  I2 <- sapply(th, function(t) { P <- plogis(dd + a * t); sum(a^2 * P * (1 - P)) })
  expect_equal(cr$info, I2, tolerance = 1e-6)
  expect_equal(cr$rel, I2 / (I2 + 1), tolerance = 1e-6)

  m <- paste("theta =~ y1 + y2 + y3 + y4", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4", "theta =~ t1 + t2 + t3 + t4", sep = "\n")
  f2 <- suppressWarnings(birt(paste(m, "V(t1 + t2 + t3 + t4) ~ theta", sep = "\n"), d, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  p <- f2$rtmb$par; S <- f2$S
  rho <- sapply(S$ld_c, function(rr) { r <- rr[S$ld$f[rr] == 1]; birt:::loading(p, S$ld[r, ]) })
  vz <- exp(2 * p$lsd2)
  IRT <- sapply(th, function(t) { s2 <- exp(2 * (p$lsig + p$kap * t)); sum(rho^2 / s2 + 2 * p$kap^2) - sum(rho / s2)^2 / (sum(1 / s2) + 1 / vz) })
  cr2 <- cond_reliability(f2, theta = th, bins = 0)
  expect_equal(cr2$info_rt, IRT, tolerance = 1e-6)
  expect_true(all(cr2$rel >= cr2$rel_acc - 1e-8) && all(cr2$se_rel > 0))
  cr3 <- cond_reliability(f2, bins = 5)
  emp <- attr(cr3, "empirical")
  expect_equal(nrow(emp), 5); expect_lt(max(abs(emp$rel_empirical - emp$rel_model)), 0.1)
  expect_output(print(cr3), "Conditional reliability")
  expect_error(cond_reliability(f2, at = c(zz = 1)), "not covariates")
})

test_that("conditional reliability: ordinal items and moderated loadings / shifts at the moderator values", {
  skip_if_not_installed("RTMB")
  ml <- rtmb_control(method = "aghq")
  hd <- sim_hgrm(N = 400, J = 5, seed = 2); hd$z <- as.numeric(scale(hd$logT))
  th <- c(-1.5, 0, 1)
  for (ty in c("grm", "gpcm", "tppcm")) {
    f <- suppressWarnings(birt(paste("f =~ y1 + y2 + y3 + y4 + y5", "f ~ z", "E(y1 + y2) ~ z", "E(y3) ~ z:f", sep = "\n"), hd,
                               ordered = TRUE, itemtype = ty, control = ml, progress = FALSE))
    for (z in c(-1, 1)) {
      cr <- cond_reliability(f, theta = th, at = c(z = z), se = FALSE, bins = 0)
      oi <- ordinal_information(f, at = list(z = z), theta = th)
      expect_equal(cr$info, oi$information, tolerance = 1e-8)                       # the item information of ordinal_information()
      nm <- sapply(th, function(t) { F <- fisher_num(f, t, xb = c(z = z)); F$info(F$all) })
      expect_equal(cr$info, nm, tolerance = 1e-6)                                     # numerical Fisher information of the likelihood
      expect_equal(cr$info_acc, cr$info)
    }
  }
  # joint model: moderated loadings and shifts of binary items and response times, a SHASH residual
  d <- sim_rtirt(N = 400, K = 5, seed = 3)
  m <- paste(rtirt_syntax(paste0("y", 1:5), paste0("t", 1:5)), "ability =~ t1 + t2", "E(y1) ~ x1", "E(y2) ~ x1:ability",
             "E(t1) ~ x1", "E(t3) ~ x1:speed", "E(t2) ~ x1:ability", sep = "\n")
  f2 <- suppressWarnings(birt(m, d, family = list(t2 = shash()), control = ml, progress = FALSE))
  cr <- cond_reliability(f2, theta = th, at = c(x1 = 1.2), se = FALSE, bins = 0)
  xb <- colMeans(f2$S$X); xb["x1"] <- 1.2
  nm <- t(sapply(th, function(t) { F <- fisher_num(f2, t, xb = xb); c(F$info(F$all), F$info(F$cat), F$info(F$rt)) }))
  expect_equal(cbind(cr$info, cr$info_acc, cr$info_rt), unname(nm), tolerance = 1e-5)
  ri <- rt_information(f2, theta = th, at = c(x1 = 1.2), limits = list(slow = .9))
  expect_equal(ri$check$info_partition + ri$check$speed_prior, cr$info_rt, tolerance = 1e-5)
})
