# SIMEX correction for MR-Egger under NOME violation

Read when `I^2_GX < 0.9` and Egger must be reported. Verbatim code block from SKILL.md "NOME violation invalidating Egger"; the trigger, mechanism and threshold stay there.

```r
library(simex)

# Precompute the weights vector rather than dividing a data-frame column in-formula:
# simex() refits the model internally on perturbed data and cannot re-evaluate
# `1 / se.outcome^2` against its own working frame, which has no se.outcome column --
# that in-formula form crashes with "object 'se.outcome' not found" inside simex()'s
# refit (simex 1.8). Precomputing the vector and passing a fully-qualified `data = dat`
# avoids it.
w <- 1 / dat$se.outcome^2
egger_lm <- lm(beta.outcome ~ beta.exposure, weights = w, data = dat, x = TRUE, y = TRUE)

egger_simex <- simex(model = egger_lm, SIMEXvariable = 'beta.exposure',
                      measurement.error = dat$se.exposure,
                      lambda = seq(0.5, 2, 0.5), B = 1000,
                      fitting.method = 'quadratic', asymptotic = FALSE)

cat('SIMEX-corrected slope:', round(coef(egger_simex)['beta.exposure'], 4), '\n')
cat('Naive Egger slope:   ', round(coef(egger_lm)['beta.exposure'], 4), '\n')
```
