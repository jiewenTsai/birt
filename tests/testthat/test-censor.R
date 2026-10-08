test_that("censor = : checks, limits on the data scale, and no censored values leave the fit unchanged", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 300, K = 4, seed = 7)
  m <- rtirt_syntax(paste0("y", 1:4), paste0("t", 1:4))
  ml <- rtmb_control(method = "aghq")
  expect_error(birt(m, d, censor = list(y1 = c(upper = 1)), control = ml, progress = FALSE), "not continuous indicators")
  expect_error(birt(m, d, censor = list(t1 = c(max = 1)), control = ml, progress = FALSE), "named number or pair")
  expect_error(birt(m, d, censor = list(t1 = c(lower = 5, upper = 4)), control = ml, progress = FALSE), "lower must be below")
  expect_error(birt(m, d, censor = list(t1 = c(upper = 60)), engine = "jags", progress = FALSE), "needs engine = \"rtmb\"")
  ds <- d; for (v in paste0("t", 1:4)) ds[[v]] <- exp(d[[v]])
  expect_error(birt(m, ds, log_rt = TRUE, censor = list(t1 = c(lower = 0)), control = ml, progress = FALSE), "positive seconds")
  # limits beyond every value: the same likelihood as without censor
  f0 <- birt(m, d, control = ml, progress = FALSE)
  expect_message(f1 <- birt(m, d, censor = list(t1 = c(upper = 100)), control = ml, progress = FALSE), "no values at or beyond")
  expect_equal(f1$rtmb$objective, f0$rtmb$objective, tolerance = 1e-10)
  # seconds with log_rt = TRUE: the limits are log-transformed with the times (same fit as on the log scale)
  up <- unname(stats::quantile(ds$t2, 0.8))
  dsc <- ds; dsc$t2 <- pmin(ds$t2, up); dlc <- d; dlc$t2 <- pmin(d$t2, log(up))
  fs <- birt(m, dsc, log_rt = TRUE, censor = list(t2 = c(upper = up)), control = ml, progress = FALSE)
  fl <- birt(m, dlc, censor = list(t2 = c(upper = log(up))), control = ml, progress = FALSE)
  expect_equal(fs$rtmb$objective, fl$rtmb$objective, tolerance = 1e-8)
  expect_equal(fs$spec$censor$n_right, sum(dsc$t2 >= up))
  expect_output(print(fs), sprintf("Censored .*t2: >= %s \\(log %s\\) %d right", format(signif(log(up), 4)), format(up), sum(dsc$t2 >= up)))
  expect_output(print(equations(fs)), "right censored")
  expect_output(print(algorithm(fs)), "censored indicators")
})

test_that("censored -logL equals brute-force integration (normal and SHASH residuals)", {
  skip_if_not_installed("RTMB")
  set.seed(4); N <- 150
  sp <- rnorm(N, 0, 0.4); e <- function() T_std(rnorm(N), 0.6, 0.8)
  d <- data.frame(t1 = 4 - sp + 0.5 * e(), t2 = 3.5 - sp + 0.4 * e(), t3 = 3.8 - sp + 0.6 * e())
  lim <- list(c(-Inf, 4.4), c(3, 4), c(3.2, Inf))
  cz <- list(t1 = c(upper = 4.4), t2 = c(lower = 3, upper = 4), t3 = c(lower = 3.2))
  for (j in 1:3) d[[j]] <- pmin(pmax(d[[j]], lim[[j]][1]), lim[[j]][2])
  for (sh in c(FALSE, TRUE)) {
    fam <- if (sh) list(t1 = shash(), t2 = shash(), t3 = shash())
    f <- suppressWarnings(birt("speed =~ -1*t1 + -1*t2 + -1*t3", d, censor = cz, family = fam, control = rtmb_control(method = "aghq"), progress = FALSE))
    p <- f$rtmb$par; s <- exp(p$lsd1)
    # independent inverse of the residual CDF (root finding) and density
    Tinv <- function(x, j) if (!sh) x else stats::uniroot(function(z) T_std(z, p$eps[j], exp(p$ldl[j])) - x, c(-60, 60), extendInt = "yes", tol = 1e-13)$root
    ldens <- function(x, j) if (!sh) dnorm(x, log = TRUE) else ld_std(x, p$eps[j], exp(p$ldl[j]))
    ll <- 0
    for (i in 1:N) {
      g <- function(sv) vapply(sv, function(v) {
        out <- dnorm(v / s, log = TRUE) - log(s)
        for (j in 1:3) { mu <- p$xi[j] - v; sdj <- exp(p$lsig[j]); y <- d[i, j]
          out <- out + if (y >= lim[[j]][2]) pnorm(Tinv((lim[[j]][2] - mu) / sdj, j), lower.tail = FALSE, log.p = TRUE)
                       else if (y <= lim[[j]][1]) pnorm(Tinv((lim[[j]][1] - mu) / sdj, j), log.p = TRUE)
                       else ldens((y - mu) / sdj, j) - log(sdj) }
        exp(out) }, 0)
      ll <- ll + log(stats::integrate(g, -8 * s, 8 * s, rel.tol = 1e-12, subdivisions = 1000)$value)
    }
    # the likelihood code (ind_ll) on a fine grid of the person latent, and the AGHQ of the fit with 81 nodes
    S <- f$S; zz <- seq(-8, 8, length.out = 4001); lg <- 0
    for (i in 1:N) { rows <- rep(i, length(zz)); v <- ind_ll(p, S, lat_values(p, S, zz, 0, rows), rows) + dnorm(zz, log = TRUE)
      lg <- lg + max(v) + log(sum(exp(v - max(v))) * diff(zz)[1]) }
    expect_equal(lg, ll, tolerance = 1e-8)
    obj <- RTMB::MakeADFun(make_aghq(S, node_data(S, centers(S, p, f$rtmb$centers$z), 81)), p, silent = TRUE)
    expect_lt(abs(obj$fn(obj$par) + ll), 1e-4)
    if (!sh) expect_lt(abs(f$rtmb$objective + ll), 1e-8)
  }
})

test_that("censored information: cond_reliability() and rt_information() match the curvature of the censored likelihood", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 500, K = 5, seed = 21)
  m <- paste("ability =~ y1 + y2 + y3 + y4 + y5", "speed =~ -1*t1 + -1*t2 + -1*t3 + -1*t4 + -1*t5", "ability =~ t1 + t2", "V(t1) ~ ability", sep = "\n")
  fam <- list(t1 = shash(), t2 = shash()); ml <- rtmb_control(method = "aghq")
  cz <- list(t1 = c(upper = 4.3), t2 = c(lower = 3.2, upper = 4.4), t3 = c(upper = 4.2))
  dc <- d; dc$t1 <- pmin(d$t1, 4.3); dc$t2 <- pmin(pmax(d$t2, 3.2), 4.4); dc$t3 <- pmin(d$t3, 4.2)
  fc <- suppressWarnings(birt(m, dc, censor = cz, family = fam, control = ml, progress = FALSE))
  th <- c(-1.5, 1)
  cr <- cond_reliability(fc, theta = th, se = FALSE, bins = 0)
  nm <- t(sapply(th, function(t) { F <- fisher_num(fc, t); c(F$info(F$all), F$info(F$rt)) }))
  expect_equal(cr$info, nm[, 1], tolerance = 1e-5)
  expect_equal(cr$info_rt, nm[, 2], tolerance = 1e-5)
  # rt_information(): limits at the item quantile u, against the censored likelihood with those limits
  f0 <- suppressWarnings(birt(m, d, family = fam, control = ml, progress = FALSE))
  ri <- rt_information(f0, theta = th, limits = list(slow = .9, fast = .05))
  p0 <- f0$rtmb$par; S0 <- f0$S
  for (side in c("slow", "fast")) {
    u <- if (side == "slow") .9 else .05
    v <- sapply(th, function(t) {
      st <- cr_setup(p0, S0, t, colMeans(S0$X))
      lim <- lapply(seq_along(S0$C), function(j) { s <- match(j, S0$sh)
        cc <- st$P0$mu[, j] + exp(st$P0$lsc[, j]) * (if (is.na(s)) qnorm(u) else T_std(qnorm(u), st$P0$e[, j], exp(p0$ldl[s])))
        if (side == "slow") c(-Inf, cc) else c(cc, Inf) })
      names(lim) <- S0$C; F <- fisher_num(f0, t, limits = lim); F$info(F$rt) })
    expect_equal(ri$limits$info_rt[ri$limits$side == side], v, tolerance = 1e-5)
  }
})

test_that("censored fits: ELGM, scores and the generated RTMB script", {
  skip_if_not_installed("RTMB")
  d <- sim_rtirt(N = 300, K = 4, seed = 3)
  m <- rtirt_syntax(paste0("y", 1:4), paste0("t", 1:4))
  dc <- d; dc$t1 <- pmin(d$t1, 4.2); dc$t2 <- pmax(d$t2, 3.4)
  cz <- list(t1 = c(upper = 4.2), t2 = c(lower = 3.4))
  f <- suppressWarnings(birt(m, dc, censor = cz, family = list(t1 = shash()), control = rtmb_control(method = "aghq"), progress = FALSE))
  e <- new.env(); e$data <- dc
  utils::capture.output(eval(parse(text = rtmb_code(f)), e))
  expect_equal(e$optA$objective, f$rtmb$objective, tolerance = 1e-5)
  S <- f$S; p <- f$rtmb$par; z1 <- rnorm(300); z2 <- rnorm(300)
  expect_equal(sum(e$joint(p, z1, z2, 1:300)), sum(ind_ll(p, S, lat_values(p, S, z1, z2, 1:300), 1:300) + dnorm(z1, log = TRUE) + dnorm(z2, log = TRUE)), tolerance = 1e-10)
  expect_lt(max(abs(colSums(estfun.birt_rtmb(f)))), 0.05)                    # casewise scores of the censored likelihood sum to ~0
  fe <- suppressWarnings(birt(m, dc, censor = cz, progress = FALSE))
  expect_true(is.finite(fe$rtmb$elgm$lognc))
  expect_output(print(fe), "Censored")
})
