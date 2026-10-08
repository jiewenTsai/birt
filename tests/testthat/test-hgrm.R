d <- sim_hgrm(N = 150, J = 5)
m <- "theta =~ y1 + y2 + y3 + y4 + y5"

test_that("moderation statements: terms, labels, anchors and identification", {
  s <- birt:::build_spec(paste(m, 'E(y1 + y2) ~ prior("ssp")*logT', "E(y3) ~ 0*logT",
                               'E(y4) ~ prior("dnorm(0, 4)")*logT       # threshold shift',
                               "E(y1 + y2 + y3 + y4 + y5) ~ pa*logT:theta", "E(theta) ~ logT", "V(theta) ~ logT", sep = "\n"), d, ordered = TRUE)
  b <- s$mod[s$mod$kind == "beta", ]; a <- s$mod[s$mod$kind == "alpha", ]
  expect_equal(b$type, c("ssp", "ssp", "none", "free")); expect_equal(b$target, c("y1", "y2", "y3", "y4"))
  expect_equal(unique(a$type), "common"); expect_equal(unique(a$label), "pa")
  expect_true("logT" %in% s$reg$x[s$reg$lhs == "theta"]); expect_equal(s$family$factor$theta$scale, "logT")
  s$dp <- dpriors(); code <- birt:::build_code(s)$code
  expect_match(code, "dlt_y4_logT ~ dnorm(0,4)", fixed = TRUE)
  expect_match(code, "lm_theta_y5[j] <- lam_theta_y5 * exp(almc_pa * logT[j])", fixed = TRUE)
  expect_match(code, "dordered.logit", fixed = TRUE)
  # hgrm() writes the same model in birt() syntax
  hs <- birt:::hgrm_syntax(m, "logT", "common", "free", c("mean", "var"), "y5")
  expect_match(hs, "E(y1 + y2 + y3 + y4) ~ logT", fixed = TRUE)
  expect_match(hs, "E(y1 + y2 + y3 + y4 + y5) ~ ca_logT*logT:theta", fixed = TRUE)   # common a: every item, also anchors
  expect_match(birt:::hgrm_syntax(m, "logT", "free", "none", "mean", "y5"), "E(y1 + y2 + y3 + y4) ~ logT:theta", fixed = TRUE)
  expect_error(birt:::hgrm_syntax(m, "logT", "ssp", "none", "mean", NULL), "must be one of")
  expect_match(hs, "V(theta) ~ logT", fixed = TRUE)
  # errors
  expect_error(birt:::build_spec(paste(m, "E(y1 + y2 + y3 + y4 + y5) ~ logT", "E(theta) ~ logT", sep = "\n"), d, ordered = TRUE), "not identified")
  expect_silent(birt:::build_spec(paste(m, "E(y1 + y2 + y3 + y4) ~ logT", "E(theta) ~ logT", sep = "\n"), d, ordered = TRUE))
  expect_error(birt:::build_spec(paste(m, "y1 ~ logT", sep = "\n"), d, ordered = TRUE), "regression outcomes must be latent")
  expect_error(birt:::build_spec(paste(m, "E(y9) ~ logT", sep = "\n"), d, ordered = TRUE), "not a latent variable or an indicator")
  expect_error(birt:::build_spec(paste(m, "V(y1) ~ logT", sep = "\n"), d, ordered = TRUE), "no residual SD")
  expect_error(birt:::build_spec(paste(m, "E(theta) ~ logT:theta", sep = "\n"), d, ordered = TRUE), "not interactions")
  expect_error(birt:::build_spec(paste(m, "V(theta) ~ logT", "V(theta) ~ logT", sep = "\n"), d, ordered = TRUE), "Duplicate|twice")
  expect_error(normal(scale = "logT"), "unused argument"); expect_error(shash(scale = "logT"), "unused argument")
  expect_error(birt:::build_spec(paste(m, "E(y1) ~ logT:Time", sep = "\n"), d, ordered = TRUE), "moderator")
  expect_error(birt:::build_spec(paste(m, "E(y1) ~ nothere", sep = "\n"), d, ordered = TRUE), "moderators not in")
  expect_error(birt:::build_spec(paste(m, "E(y1) ~ 0.5*logT", sep = "\n"), d, ordered = TRUE), "fixed nonzero")
  expect_error(birt:::build_spec("f1 =~ y1 + y2 + y3\nf2 =~ y4\nf1 =~ y5", d, ordered = TRUE), "first line")
  expect_error(birt:::build_spec("f1 =~ y1 + y2 + y3\nf2 =~ y4 + y5\nf1 =~ y4", d, ordered = TRUE), "exactly one latent variable")
  expect_error(birt:::build_spec(m, d, ordered = TRUE, family = list(y1 = normal())), "itemtype")
  expect_error(birt:::hgrm_syntax("f =~ y1 + y2\ng =~ y3 + y4", "logT", "none", "none", character(0), NULL), "one-factor")
})

test_that("ordinal indicators with response times and binary items in one model", {
  skip_if_not_installed("RTMB")
  dd <- d; set.seed(3); dd$t1 <- -0.3 * scale(dd$logT)[, 1] + rnorm(150, 0, 0.5); dd$t2 <- rnorm(150); dd$x1 <- rbinom(150, 1, 0.5)
  s <- birt:::build_spec("theta =~ y1 + y2 + y3 + x1\nspeed =~ t1 + t2\nE(t1) ~ logT\nV(t2) ~ logT", dd, ordered = c("y1", "y2", "y3"))
  expect_equal(s$ordinal, c("y1", "y2", "y3")); expect_equal(s$binary, "x1"); expect_equal(s$cont, c("t1", "t2"))
  expect_equal(s$family$ind$t2$scale, "logT")
  f <- suppressWarnings(birt("theta =~ y1 + y2 + y3 + x1\nspeed =~ t1 + t2\nE(t1) ~ logT\nV(t2) ~ logT", dd, ordered = c("y1", "y2", "y3"),
                             engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  expect_equal(f$rtmb$convergence, 0)
  e <- f$estimates
  expect_true(all(c("E(t1)", "V(t2)") %in% e$lhs)); expect_equal(sum(e$op == "|"), 9)
  # the marginal likelihood equals brute-force integration over (theta, speed)
  S <- f$S; p <- f$rtmb$par
  g <- seq(-7, 7, length.out = 161); h <- diff(g[1:2]); G <- expand.grid(z1 = g, z2 = g); total <- 0
  for (i in seq_len(S$N)) {
    rows <- rep(i, nrow(G))
    ll <- birt:::ind_ll(p, S, birt:::lat_values(p, S, G$z1, G$z2, rows), rows) + dnorm(G$z1, log = TRUE) + dnorm(G$z2, log = TRUE)
    total <- total + max(ll) + log(sum(exp(ll - max(ll))) * h^2)
  }
  expect_equal(-f$rtmb$objective, total, tolerance = 1e-4)
})

test_that("a moderated GRM fits by JAGS; the marginal likelihood matches integration", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  fit <- birt(paste(m, 'E(y1 + y2 + y3 + y4 + y5) ~ prior("ssp")*logT', "E(y1 + y2 + y3 + y4 + y5) ~ ca_logT*logT:theta",
                    "E(theta) ~ logT", "V(theta) ~ logT", sep = "\n"), d, ordered = TRUE,
              n_chains = 2, n_iter = 400, n_burn = 200, progress = FALSE, engine = "jags")
  expect_s3_class(fit, "birt")
  expect_equal(nrow(scores(fit)), 150); expect_true(!is.null(fit$loo))
  expect_equal(nrow(fit$dif), 10); expect_true(all(c("p_incl", "BF10") %in% names(fit$dif)))
  expect_equal(nrow(fit$impact), 2); expect_equal(nrow(fit$precision), 3)
  P <- birt:::prep_pars(fit, 3); LL <- birt:::loglik_compute(P, birt:::obs_matrices(fit))
  s <- 2; i <- 4; M <- P$M; z <- fit$data$logT[i]; sp <- fit$spec
  mu <- M[s, "beta_theta_logT"] * z; sdv <- exp(M[s, "psi_theta_logT"] * z); al <- M[s, "almc_ca_logT"]
  fth <- function(th) sapply(th, function(t0) { l <- 0
    for (v in sp$ordinal) { a <- M[s, sprintf("lam_theta_%s", v)] * exp(al * z); dl <- M[s, sprintf("dlt_%s_logT", v)] * z
      b <- M[s, sprintf("b_%s[%d]", v, seq_len(sp$ord$K[match(v, sp$ord$items)] - 1))]
      up <- c(1, plogis(a * (t0 - dl - b)), 0); y <- sp$ord$Y[i, v]; l <- l + log(up[y] - up[y + 1]) }
    exp(l) * dnorm(t0, mu, sdv) })
  expect_equal(LL[s, i], log(integrate(fth, -12, 12, rel.tol = 1e-12)$value), tolerance = 1e-6)
  u <- as.character(equations(fit)); expect_true(any(grepl("Priors", u))); expect_true(any(grepl("spike-and-slab", u)))
  pc <- ppc(fit, 10); expect_s3_class(pc, "birt_ppc"); expect_true(is.data.frame(pc$categories))
  expect_true(is.data.frame(ordinal_information(fit)))
  sm <- summary(fit); expect_output(print(sm), "Moderation"); expect_s3_class(sm$moderation, "data.frame")
  expect_true(all(sm$selection$evidence %in% c("strong for 0", "moderate for 0", "inconclusive", "moderate for nonzero", "strong for nonzero")))
  pi_hat <- mean(as.matrix(fit$draws)[, "p_incl_dlt_logT"])          # learned inclusion probability of the shifts
  sel <- sm$selection; expect_equal(sel$prior_incl, rep(pi_hat, 5))
  expect_equal(sel$BF10, (sel$p_incl / (1 - sel$p_incl)) / (pi_hat / (1 - pi_hat)))
  g0 <- birt(m, d, ordered = TRUE, n_chains = 2, n_iter = 300, n_burn = 150, progress = FALSE, engine = "jags")
  expect_null(g0$dif); expect_true(is.data.frame(ordinal_information(g0))); expect_s3_class(ppc(g0, 5), "birt_ppc")
  s0 <- summary(g0); expect_output(print(s0), "[$]items"); expect_true("t3" %in% names(s0$items))
})

test_that("ordinal models by maximum likelihood (RTMB): direct integration and the earlier pipeline", {
  skip_if_not_installed("RTMB")
  mm <- paste(m, "E(y1 + y2 + y3 + y4 + y5) ~ pa*logT:theta", "E(y1) ~ logT", "E(theta) ~ logT", "V(theta) ~ logT", sep = "\n")
  r <- suppressMessages(birt(mm, d, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  S <- r$S; p <- r$rtmb$par
  expect_equal(-r$rtmb$objective, sum(vapply(seq_len(S$N), function(i) {
    z <- S$X[i, "logT"]; mu <- p$beta1 * z; sdv <- exp(p$alpha1 * z)
    log(integrate(function(t) sapply(t, function(th) exp(birt:::ind_ll(p, S, list(th), i)) * dnorm(th, mu, sdv)),
                  mu - 12 * sdv, mu + 12 * sdv, rel.tol = 1e-10)$value) }, 0)), tolerance = 1e-6)
  expect_true(all(!is.na(r$estimates$se[!r$estimates$fixed])))
  expect_equal(nrow(r$dif), 6); expect_true(all(c("z", "pvalue") %in% names(r$dif)))
  expect_output(print(summary(r)), "Moderation \\(E: expectation")
  # the same model and data as the earlier ordinal pipeline (scripts/80_birt_baseline.R): -log L 1441.958001
  d3 <- sim_hgrm(N = 300, J = 5)
  r3 <- suppressMessages(birt(mm, d3, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  expect_equal(r3$rtmb$objective, 1441.958001, tolerance = 1e-5 / 1441)
  g <- birt(m, d, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
  expect_output(compare(GRM = g, mod = r), "dAIC")
  expect_error(birt(paste(m, 'E(y1) ~ prior("ssp")*logT', sep = "\n"), d, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "jags")
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_true(is.data.frame(ordinal_information(r)))
})

test_that("ordinal RTMB: two moderators in equations(), convergence(); ELGM", {
  skip_if_not_installed("RTMB")
  dd <- d; dd$age <- stats::rnorm(nrow(dd))
  r <- suppressMessages(birt(paste(m, "E(y1 + y2 + y3 + y4 + y5) ~ pa*logT:theta", "E(y1) ~ age", "E(theta) ~ logT + age", "V(theta) ~ age", sep = "\n"),
                             dd, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  expect_true(any(grepl("age", as.character(equations(r)))))
  expect_output(convergence(r), "nlminb_code")
  e <- suppressWarnings(birt(paste(m, "E(theta) ~ logT", sep = "\n"), d, ordered = TRUE, engine = "rtmb",
                             control = rtmb_control(method = "elgm"), progress = FALSE))
  expect_true(all(c("y1|t1", "theta~logT") %in% paste0(e$estimates$lhs, e$estimates$op, e$estimates$rhs)))
  expect_true(all(nzchar(e$estimates$prior[!e$estimates$fixed & !e$estimates$op %in% c("irt_a", "irt_b")])))
})

test_that("partial credit models: category probabilities, ML against direct integration and mirt", {
  th <- seq(-2, 2, 0.5); b <- c(-1, 0.3, 1.2)
  for (ty in c("grm", "gpcm", "tppcm")) expect_equal(rowSums(birt:::cat_probs(ty, 1.3, th, b, c(1, 0.7, 1.4))), rep(1, length(th)))
  expect_equal(birt:::cat_probs("gpcm", 1.3, th, b), birt:::cat_probs("tppcm", 1.3, th, b, c(1, 1, 1)))
  expect_equal(birt:::cat_probs("gpcm", 1.3, th, 0.4)[, 2], plogis(1.3 * (th - 0.4)))               # two categories: 2PL
  skip_if_not_installed("RTMB")
  set.seed(5); N <- 400; J <- 5; th <- rnorm(N)
  a <- runif(J, 0.8, 2); b <- t(apply(matrix(rnorm(J * 3), J), 1, sort)); r <- cbind(1, matrix(runif(J * 2, 0.5, 1.6), J))
  Y <- sapply(1:J, function(j) { P <- birt:::cat_probs("tppcm", a[j], th, b[j, ], r[j, ]); 1 + rowSums(runif(N) > t(apply(P, 1, cumsum))[, -4]) })
  dd <- data.frame(Y); names(dd) <- paste0("y", 1:J); mo <- paste("f =~", paste(names(dd), collapse = " + "))
  for (ty in c("gpcm", "tppcm")) {
    f <- birt(mo, dd, ordered = TRUE, itemtype = ty, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
    S <- f$S; p <- f$rtmb$par
    bf <- sum(vapply(seq_len(N), function(i) log(integrate(function(t) sapply(t, function(x) exp(birt:::ind_ll(p, S, list(x), i)) * dnorm(x)),
                                                           -12, 12, rel.tol = 1e-10)$value), 0))
    expect_equal(-f$rtmb$objective, bf, tolerance = 1e-6)
    expect_equal(sum(f$estimates$op == "|a"), if (ty == "tppcm") 2 * J else 0)
    out <- capture.output(print(summary(f)$parameters))           # section notes only where they apply
    expect_true(any(grepl(sprintf("^Thresholds \\(partial credit items: step difficulties%s\\):", if (ty == "tppcm") "; [|]a: step discriminations" else "\\)?"), out)))
    expect_false(any(grepl("SHASH", out)))
  }
  skip_if_not_installed("mirt")
  g <- birt(mo, dd, ordered = TRUE, itemtype = "gpcm", engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))
  mg <- mirt::mirt(dd, 1, itemtype = "gpcm", verbose = FALSE, TOL = 1e-7)
  expect_equal(2 * g$rtmb$objective, -2 * mg@Fit$logLik, tolerance = 1e-3 / 1e4)
  cf <- mirt::coef(mg, IRTpars = TRUE, simplify = TRUE)$items; e <- g$estimates
  expect_lt(max(abs(e$est[e$op == "=~"] - cf[, "a"])), 1e-3)
  expect_lt(max(abs(matrix(e$est[e$op == "|"], J, byrow = TRUE) - cf[, c("b1", "b2", "b3")])), 1e-3)
  tp <- birt(mo, dd, ordered = TRUE, itemtype = "tppcm", engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))   # the saturated TPPCM is a nominal model
  mn <- mirt::mirt(dd, 1, itemtype = "nominal", verbose = FALSE, TOL = 1e-8, technical = list(NCYCLES = 5000))
  expect_equal(2 * tp$rtmb$objective, -2 * mn@Fit$logLik, tolerance = 1e-2 / 1e4)
  expect_error(birt(mo, dd, ordered = TRUE, itemtype = "pcm", engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")), "itemtype must be one of")
})

test_that("partial credit models by JAGS: the marginal likelihood matches integration", {
  skip_on_cran()
  skip_if(Sys.which("jags") == "", "JAGS not installed")
  set.seed(6); N <- 200; th <- rnorm(N); z <- rnorm(N)
  Y <- sapply(1:4, function(j) { P <- birt:::cat_probs("tppcm", 1.2 * exp(0.2 * z), th, c(-1, 0, 1), c(1, 1.3, 0.8)); 1 + rowSums(runif(N) > t(apply(P, 1, cumsum))[, -4]) })
  dd <- data.frame(Y, z = z); names(dd)[1:4] <- paste0("y", 1:4)
  mo <- "f =~ y1 + y2 + y3 + y4\nE(y1 + y2 + y3 + y4) ~ pa*z:f\nE(y1) ~ z"
  for (ty in c("gpcm", "tppcm")) {
    f <- birt(mo, dd, ordered = paste0("y", 1:4), itemtype = ty, n_chains = 2, n_iter = 400, n_burn = 200, progress = FALSE, engine = "jags")
    P <- birt:::prep_pars(f, 2); LL <- birt:::loglik_compute(P, birt:::obs_matrices(f)); M <- P$M; s <- 1; i <- 5; zi <- dd$z[i]
    fth <- function(t) sapply(t, function(t0) { l <- 0
      for (v in paste0("y", 1:4)) { a <- M[s, sprintf("lam_f_%s", v)] * exp(M[s, "almc_pa"] * zi); dl <- if (v == "y1") M[s, "dlt_y1_z"] * zi else 0
        rr <- if (ty == "tppcm") c(1, M[s, sprintf("r_%s[%d]", v, 2:3)]) else rep(1, 3)
        l <- l + log(birt:::cat_probs(ty, a, t0 - dl, M[s, sprintf("b_%s[%d]", v, 1:3)], rr)[1, f$spec$ord$Y[i, v]]) }
      exp(l) * dnorm(t0) })
    expect_equal(LL[s, i], log(integrate(fth, -12, 12, rel.tol = 1e-12)$value), tolerance = 1e-6)
    expect_true(is.data.frame(ordinal_information(f))); expect_s3_class(ppc(f, 5), "birt_ppc")
    expect_true(any(grepl(if (ty == "gpcm") "generalized partial credit" else "two-parameter partial credit", as.character(equations(f)))))
  }
})

test_that("hgrm() fits by RTMB: defaults, DIF with anchors, same fit as its syntax", {
  skip_if_not_installed("RTMB")
  h <- hgrm(m, d, moderators = "logT", progress = FALSE, control = rtmb_control(method = "aghq"))                       # a = "common", impact = "mean"
  expect_s3_class(h, "birt_rtmb"); expect_equal(h$likelihood, "marginal")
  expect_match(h$model, "E(y1 + y2 + y3 + y4 + y5) ~ ca_logT*logT:theta", fixed = TRUE)
  h2 <- hgrm(m, d, moderators = "logT", a = "none", b = "free", anchor = c("y4", "y5"), progress = FALSE, control = rtmb_control(method = "aghq"))
  expect_equal(nrow(h2$dif), 3); expect_true(all(c("z", "pvalue") %in% names(h2$dif)))
  dz <- d; dz$logT <- as.numeric(scale(dz$logT))
  b2 <- birt(h2$model, dz, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq"))  # hgrm() = its syntax on standardized moderators
  expect_equal(b2$rtmb$objective, h2$rtmb$objective, tolerance = 1e-8)
})

test_that("summaries of a GRM do not print partial credit or SHASH notes", {
  skip_if_not_installed("RTMB")
  dz <- d; dz$z <- as.numeric(scale(dz$logT))
  f <- suppressWarnings(birt(paste(m, "V(theta) ~ z", sep = "\n"), dz, ordered = TRUE, progress = FALSE, control = rtmb_control(method = "aghq")))
  out <- capture.output(print(summary(f)$parameters))
  expect_true("Thresholds:" %in% out); expect_true("Distributions:" %in% out)
  expect_false(any(grepl("partial credit|SHASH", out)))
})

test_that("BF10 evidence categories include BF10 = Inf", {
  expect_equal(birt:::bf_evidence(c(0, 0.05, 0.2, 1, 5, 10, Inf, NA)),
               c("strong for 0", "strong for 0", "moderate for 0", "inconclusive", "moderate for nonzero", "strong for nonzero", "strong for nonzero", NA))
})

test_that("RTMB: V() effects are unbounded and ordinal loadings are in the measurement block", {
  skip_if_not_installed("RTMB")
  m <- "theta =~ y1 + y2 + y3 + y4 + y5"
  d <- sim_hgrm(N = 300, J = 5, seed = 4); d$z <- as.numeric(scale(d$logT)); d$z20 <- d$z / 20
  r1 <- suppressMessages(birt(paste(m, "V(theta) ~ z", sep = "\n"), d, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  r2 <- suppressMessages(birt(paste(m, "V(theta) ~ z20", sep = "\n"), d, ordered = TRUE, engine = "rtmb", progress = FALSE, control = rtmb_control(method = "aghq")))
  expect_equal(r1$rtmb$objective, r2$rtmb$objective, tolerance = 1e-6)   # rescaling the moderator: same likelihood
  sc <- estfun.birt_rtmb(r1)
  pa <- attr(sc, "parameters"); expect_true(all(pa$block[grepl("^l(pos|real)\\[", pa$par)] == "measurement"))
})
