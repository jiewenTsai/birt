# Example 6. Likert scale + one total response time: moderated graded response model (MNLFA)
#
# E(): the conditional expectation (linear predictor) of an item or of the trait; V(): the
# variance of the trait. E(item) ~ z shifts the thresholds (uniform DIF), E(item) ~ z:SA
# moderates the discrimination (nonuniform DIF). The models are fitted by ELGM (the default);
# the score-based DIF screen needs maximum likelihood, control = rtmb_control(method = "aghq").
library(birt)
d <- sim_hgrm(N = 800, J = 8)                    # y1..y8 (1-4), Time, logT; true DIF (threshold shift) in y1, y2
d$z <- as.numeric(scale(d$logT))                 # effects per SD (birt() uses moderators as given)
items <- "SA =~ y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8"
all8 <- "y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8"
precision <- paste(items, sprintf("E(%s) ~ pa*z:SA", all8), "E(SA) ~ z", sep = "\n")   # one common precision effect (equal labels)

# 1. ELGM: graded response model and the precision model
grm  <- birt(items, d, ordered = TRUE)
prec <- birt(precision, d, ordered = TRUE)
compare(GRM = grm, precision = prec)             # log marginal likelihoods, Bayes factors

# 2. which items' parameters change with z? one maximum likelihood fit, every item (Holm-adjusted)
prec_ml <- birt(precision, d, ordered = TRUE, control = rtmb_control(method = "aghq"))
score_test(prec_ml, d$z, by_item = TRUE)         # y1 and y2 have true DIF; other flags are screening hints

# 3. ELGM again: free threshold shifts for the flagged items (the others are anchors)
full <- birt(paste(precision, "E(y1 + y2) ~ z", "V(SA) ~ z", sep = "\n"), d, ordered = TRUE)
compare(GRM = grm, precision = prec, precision_DIF = full)
summary(full)                     # item parameters, moderation effects with 95% intervals, impact
summary(full)$precision           # precision of the scores by tercile of z: fast respondents measured less precisely
equations(full)
oi <- ordinal_information(full)   # test information at the 10th / 50th / 90th percentile of z
oi[oi$theta %in% c(-2, 0, 2), ]
plot(full, type = "information")
plot(full, type = "dif")          # expected item scores at the 10th / 90th percentile of z

# 4. hgrm(): the same models from a few options (it standardizes the moderators; the syntax it writes is fit$model)
h <- hgrm(items, d, moderators = "logT", a = "common", b = "free", anchor = paste0("y", 3:8), impact = c("mean", "var"), label = "DIF")
cat(h$model)

# 5. nonlinear moderation: a piecewise term with a knot at 2 seconds per item (log scale)
d$logT_fast <- pmin(log(d$Time / 196) - log(2), 0)   # 196 items in the whole questionnaire
hp <- hgrm(items, d, moderators = "logT", label = "precision")   # defaults: a = "common", b = "none", impact = "mean"
pw <- hgrm(items, d, moderators = c("logT", "logT_fast"), a = c("common", "common"), impact = "mean", label = "piecewise")
compare(precision = hp, piecewise = pw)

# 6. partial credit instead of graded response (experimental with moderation)
gp <- birt(paste(precision, "E(y1 + y2) ~ z", sep = "\n"), d, ordered = TRUE, itemtype = "gpcm")
compare(GRM = full, GPCM = gp)

# 7. Bayesian DIF screening with spike-and-slab priors on every shift and discrimination (JAGS)
if (requireNamespace("rjags", quietly = TRUE)) {
  ssp <- birt(paste(items,
                    sprintf('E(%s) ~ prior("ssp")*z', all8),        # threshold shifts (uniform DIF)
                    sprintf('E(%s) ~ prior("ssp")*z:SA', all8),     # discriminations (nonuniform DIF)
                    "E(SA) ~ z", "V(SA) ~ z", sep = "\n"), d, ordered = TRUE, engine = "jags", progress = FALSE)
  print(summary(ssp)$selection)   # p_incl, prior_incl (learned inclusion probability) and BF10 per effect
  convergence(ssp)                # these models mix slowly: check R-hat, use long chains
}
