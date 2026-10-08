test_that("dif_tree() keeps a single node without DIF and splits on a DIF variable", {
  skip_on_cran(); skip_if_not_installed("partykit"); skip_if_not_installed("RTMB"); skip_if_not_installed("strucchange")
  m <- rtirt_syntax(paste0("y", 1:6), paste0("t", 1:6))
  d <- sim_rtirt(N = 600, K = 6, seed = 4); set.seed(4)
  t0 <- dif_tree(m, d, data.frame(g = factor(sample(1:2, 600, TRUE))))
  expect_equal(length(partykit::nodeids(t0)), 1)
  g <- rep(1:2, each = 300); d1 <- d
  d1$y1[g == 2] <- rbinom(300, 1, 0.92); d1$y2[g == 2] <- rbinom(300, 1, 0.08)       # strong DIF in y1, y2 for group 2
  t1 <- dif_tree(m, d1, data.frame(g = factor(g), noise = rnorm(600)))
  expect_gt(length(partykit::nodeids(t1)), 1)
  expect_equal(as.character(partykit::split_node(partykit::node_party(t1))$varid |> (\(i) names(t1$data)[i])()), "g")
})

test_that("dif_tree() cuts numeric partitioning variables into quantile groups and finds a numeric DIF split", {
  g <- birt:::quantile_groups(c(1:100, 100, 100), 10)
  expect_true(is.ordered(g)); expect_lte(nlevels(g), 10); expect_true(all(table(g) >= 9))
  skip_on_cran(); skip_if_not_installed("partykit"); skip_if_not_installed("RTMB"); skip_if_not_installed("strucchange")
  set.seed(4); N <- 600; J <- 6
  th <- rnorm(N); x <- runif(N, 0, 10); b <- seq(-1, 1, length.out = J)
  B <- matrix(b, N, J, byrow = TRUE); B[x > 6, 1] <- B[x > 6, 1] + 1; B[x > 6, 2] <- B[x > 6, 2] - 1   # DIF in y1, y2 for x > 6
  Y <- (matrix(runif(N * J), N, J) < plogis(1.2 * (th - B))) * 1; colnames(Y) <- paste0("y", 1:J)
  tr <- dif_tree(paste("theta =~", paste(colnames(Y), collapse = " + ")), as.data.frame(Y), data.frame(x = x, noise = rnorm(N)))
  sp <- partykit::split_node(partykit::node_party(tr))
  expect_equal(names(tr$data)[sp$varid], "x")
  cut_at <- as.numeric(sub(".*,([0-9.]+)\\]$", "\\1", levels(tr$data$x)[sp$breaks]))   # upper end of the last group on the left
  expect_lt(abs(cut_at - 6), 0.75)
  expect_error(dif_tree("theta =~ y1 + y2", as.data.frame(Y), data.frame(x = x), nbins = 1), "nbins")
})
