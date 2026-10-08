# Example 9. The van der Linden correlation is a common cross-relation (Cholesky)
#
# speed = c * ability + e (c = r * sd(speed), Var(e) = (1 - r^2) Var(speed)) turns
#   log t_ij = xi_i - speed_j + eps_ij,  Corr(ability, speed) = r
# into
#   log t_ij = xi_i - e_j - c * ability_j + eps_ij,  e independent of ability,
# a cross-loading of ability that is the same for every item i. The cross-relation model of
# Bolsinova and Tijmstra lets it differ by item; their mean carries the person-level
# correlation, the deviations the item-level conditional dependence.
library(birt)
d <- sim_rtirt(N = 800, K = 6, seed = 3)
it <- paste0("y", 1:6); tm <- paste0("t", 1:6)
base <- paste0("ability =~ ", paste(it, collapse = " + "), "\nspeed =~ ", paste(sprintf("-1*%s", tm), collapse = " + "))
m_hier  <- paste(base, "ability ~~ speed", sep = "\n")                                          # correlation
m_equal <- paste(base, paste0("ability =~ ", paste(sprintf("r*%s", tm), collapse = " + ")), sep = "\n")  # equal cross-loadings
m_cross <- paste(base, paste0("ability =~ ", paste(tm, collapse = " + ")), sep = "\n")          # item-specific

# maximum likelihood: the first two models have the same likelihood
ml <- lapply(list(correlation = m_hier, equal_cross = m_equal, item_specific = m_cross),
             function(m) birt(m, d, control = rtmb_control(method = "aghq")))
sapply(ml, logLik)                                     # correlation = equal_cross
compare(ml)
2 * (logLik(ml$item_specific) - logLik(ml$equal_cross))  # LR test of item-level dependence, df = K - 1

# ELGM: the same likelihood with different priors (r ~ U(-1, 1) vs a cross-loading ~ N(0, 1))
# gives different marginal likelihoods
el <- lapply(list(correlation = m_hier, equal_cross = m_equal, item_specific = m_cross),
             function(m) birt(m, d))
compare(el)
equations(el$equal_cross)
