#' birt: BIvariate Response Time Modeling with 'lavaan'-Style Syntax
#'
#' Main function [birt()]: by default RTMB with ELGM (approximate Bayes on the marginal likelihood);
#' maximum likelihood with `rtmb_control(method = "aghq")`; the optional JAGS engine with
#' `engine = "jags"`. Syntax helper [rtirt_syntax()];
#' residual distributions [ald()], [normal()], [shash()]; results via `summary()`,
#' [estimates()], [scores()], [reliability()], [fit_indices()],
#' [quantiles()], [ppc()], [convergence()], [compare()], `plot()`; JAGS code via [jags_code()].
#' See the tutorial (`system.file("tutorial", "birt_tutorial.qmd", package = "birt")`).
#'
#' The methods of the generics shared with birtRcpp (`estimates`, `scores`, `reliability`,
#' `fit_indices`, `convergence`) are also exported, so they dispatch whichever package was
#' attached last.
#' @keywords internal
#' @importFrom coda as.mcmc.list
#' @rawNamespace export(estimates.birt, estimates.birt_rtmb, scores.birt_fit, scores.birt_rtmb, reliability.birt_fit, reliability.birt_rtmb, fit_indices.birt_fit, fit_indices.birt_rtmb, convergence.birt_fit, convergence.birt_rtmb)
"_PACKAGE"
