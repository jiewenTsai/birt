# Example 8. The RTMB engine: ELGM and maximum likelihood on the same marginal likelihood,
# and the maximum likelihood tools (score tests, influence functions, DIF trees, RTMB code)
library(birt)
d <- sim_rtirt(N = 500, K = 10)
y <- paste0("y", 1:10); t <- paste0("t", 1:10)
model <- rtirt_syntax(y, t)                                 # ability ~~ speed

fe <- birt(model, d)                                        # ELGM (default): posterior, logML
fm <- birt(model, d, control = rtmb_control(method = "aghq"))   # maximum likelihood: estimates, SE, AIC / BIC
cbind(ELGM = coef(fe), ML = coef(fm))[1:6, ]
fit_indices(fe); fit_indices(fm)

# score-based tests of measurement invariance (maximum likelihood fits only)
score_test(fm, d$x1 > 0, by_item = TRUE, impact_check = FALSE)   # a split unrelated to the items
score_test(fm, rowSums(d[t]))                                   # along the total log time

# influence functions of quantile effects: robust SEs and invariance of one effect
m5 <- paste(rtirt_syntax(y, t, cross = "free"),
            sprintf("V(%s) ~ ability", paste(t, collapse = " + ")), sep = "\n")
f5 <- birt(m5, d, control = rtmb_control(method = "aghq"))
qs <- quantile_score(f5, p = c(0.1, 0.9))
qs                                   # est, model-based se, sandwich se.robust
score_test(qs, d$x2)                 # does an effect change along x2?

# DIF tree (needs partykit): y3 made easier for x2 > 0
if (requireNamespace("partykit", quietly = TRUE)) {
  d2 <- d; set.seed(5); g <- d2$x2 > 0
  d2$y3[g] <- pmax(d2$y3[g], rbinom(sum(g), 1, 0.3))
  z <- data.frame(x1_high = factor(d2$x1 > 0), x2_high = factor(d2$x2 > 0))
  tr <- dif_tree(paste("ability =~", paste(y, collapse = " + ")), d2, z)
  print(tr)
  # numeric partitioning variables work too, but the split search refits the model at every
  # candidate cut point and takes minutes
}

# a stand-alone RTMB script of the fitted model
invisible(capture.output(code <- rtmb_code(fm)))
cat(head(strsplit(code, "\n")[[1]], 20), sep = "\n")
