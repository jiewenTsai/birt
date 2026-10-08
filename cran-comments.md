## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

(Locally, `R CMD check --as-cran` also reports "unable to verify current time" and the
HTML Tidy version, both from the check environment.)

## Test environments

* local macOS (aarch64), R 4.5.3
* win-builder (devel and release): to be run before submission

## Notes for the reviewer

* The main engine is 'RTMB' (Imports). Its examples that fit models are in `\donttest{}`.
* The JAGS engine (`engine = "jags"`) is optional: it needs the JAGS library
  (SystemRequirements) and 'R2jags'/'rjags' (Suggests). Its examples are guarded by
  `requireNamespace()` and its tests are skipped when JAGS is not available.
* 'partykit', 'strucchange' and 'sandwich' (Suggests) are only used by `dif_tree()` and the
  score tests, which check for them.
* JAGS sampling runs in a background R session ('callr') only when `progress = TRUE`;
  examples and tests use `progress = FALSE`.
