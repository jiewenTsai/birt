# Example 2. Model comparison by log marginal likelihood (ELGM) and by AIC / BIC (maximum likelihood)
library(birt)
d <- sim_rtirt(N = 500, K = 10)        # true cross-loadings: 0 for odd items, -0.3 for even items
y <- paste0("y", 1:10); t <- paste0("t", 1:10)
acc <- paste("ability =~", paste(y, collapse = " + "))
spd <- paste("speed =~", paste0("-1*", t, collapse = " + "))
models <- list(
  independent = paste(acc, spd, sep = "\n"),                                                    # no relation
  correlation = paste(acc, spd, "ability ~~ speed", sep = "\n"),                                # van der Linden
  equal_cross = paste(acc, spd, paste("ability =~", paste0("r*", t, collapse = " + ")), sep = "\n"),  # one shared cross-loading
  item_cross  = paste(acc, spd, paste("ability =~", paste(t, collapse = " + ")), sep = "\n"))        # item cross-loadings

# ELGM (default): log marginal likelihoods and Bayes factors
fits <- lapply(names(models), function(m) birt(models[[m]], d, label = m))
names(fits) <- names(models)
compare(fits)                            # logML, dlogML, BF_best (they depend on the priors)
sapply(fits, reliability)                # reliability of ability and speed

# the same comparison with BSEM-type small-variance priors on the cross-loadings (SD 0.1)
bsem <- birt(models$item_cross, d, dp = dpriors(cross = "dnorm(0, 100)"), label = "item_cross, SD 0.1")
compare(item_cross = fits$item_cross, item_cross_bsem = bsem)

# maximum likelihood on the same marginal likelihood: AIC and BIC
ml <- lapply(models, function(m) birt(m, d, control = rtmb_control(method = "aghq")))
compare(ml)                              # -2logL, AIC, BIC and their differences

# spike-and-slab selection of the cross-loadings needs JAGS: see ex05 / the tutorial, section 9
