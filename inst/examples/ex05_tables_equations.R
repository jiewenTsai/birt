# Example 5. Summary tables as objects, the model as equations, the computation as steps
library(birt)
d <- sim_rtirt(N = 500, K = 8)
model <- rtirt_syntax(paste0("y", 1:8), paste0("t", 1:8), cross = "free")
fit <- birt(model, d)       # ELGM

s <- summary(fit)          # computed once; printed under the names of its tables
names(s)                   # header, fit, items, se, parameters, moderation, precision, reliability
s$items                    # one row per item (loadings, intercepts, residual variances, IRT a / b)
s$se                       # posterior SDs in the same layout
print(s, tables = c("fit", "parameters"))     # some tables only

# export: the long parameter table (blavaan-like columns, with the prior of each parameter)
write.csv(as.data.frame(s), file.path(tempdir(), "parameters.csv"), row.names = FALSE)
write.csv(s$items, file.path(tempdir(), "items.csv"), row.names = FALSE)

# the fitted model written out as equations: Unicode for the console, LaTeX for a paper
equations(fit)
tex <- file.path(tempdir(), "model.tex")
writeLines(as.character(equations(fit, "latex")), tex)   # \[ \input{model.tex} \] in the document

# how the ELGM fit was computed (inner AGHQ, posterior mode, outer grid, log marginal likelihood)
algorithm(fit)

# a stand-alone RTMB script of the model for an appendix (the maximum likelihood script,
# started at the posterior mode; it is printed, so capture it and write it to a file)
rfile <- file.path(tempdir(), "model_rtmb.R")
invisible(capture.output(code <- rtmb_code(fit, file = rfile)))
length(readLines(rfile))
