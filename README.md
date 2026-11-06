# rrcsurv

Risk set regression calibration (RRC) for Cox proportional hazards and linear excess relative risk
(ERR) models with an error-prone cumulative exposure, following Yang, Spiegelman and Wang.

At each event time t, the measurement error model is fitted to the validation risk set: the pairs
of the validation subjects still at risk at t whose dose is part of the cumulative exposure at t.
Every surrogate dose of the subjects at risk is replaced by its calibrated value, and the cumulative
exposure is rebuilt for that risk set. The sandwich variance accounts for the estimated calibration
parameters.

- Risk models: Cox, `exp(beta x + gamma'Z)`, and linear ERR, `exp(gamma'Z) (1 + beta x)`.
- Measurement error models: additive (`c ~ ...`) or multiplicative (`log(c) ~ ...`, lognormal),
  with any covariates on the right-hand side.
- Internal or external validation study; error-free confounders; lagged exposure.
- Variance: sandwich (the paper's, or a joint version for internal validation) or bootstrap.
- `cumfit()`: the same Cox / ERR models without correction, for an exposure measured without
  error or for a naive analysis, with model-based and robust (clustered) standard errors.

The package is in development.

## Installation

```r
# install.packages("remotes")
remotes::install_github("fyGuo/rrcsurv")
```

## Usage

The main study is in person-time format, with one row per subject and interval `(start, stop]`. The
`dose` column holds the dose received in that interval (not the running total). The `surrogate`
flag marks rows whose dose is the error-prone value; rows whose true dose is known, such as
validated years in an internal validation study, carry the true value with `surrogate = FALSE`.
The validation data have one row per (true, surrogate) pair, with the surrogate in a column of the
same name as the exposure. A subject may have pairs at several times; give each pair's dose time
(the end of the dose's interval, on the main study's time scale) in `dose_time`, so that the
calibration at event time t uses only the doses received before t (before t - lag with a lag).

```r
library(rrcsurv)
library(survival)

fit <- rrc(Surv(start, stop, event) ~ rrcx(dose, surrogate = is_estimated, lag = 5) + smoke,
           data = main, id = "id",
           validation = validation,              # columns id, c (true dose), dose (surrogate), year, ...
           dose_time = "year",                   # when each validated dose was received
           calibration = log(c) ~ log(dose),     # multiplicative error; use c ~ dose for additive
           model = "err")
summary(fit)
confint(fit)
```

The uncorrected fit takes the same formula: `cumfit(Surv(start, stop, event) ~ rrcx(dose) + smoke, data = main, id = "id", model = "err")`.

See `?rrc`, `?"rrc-methods"` and `?cumfit`.

## Validation

Without `dose_time`, the estimates and sandwich standard errors reproduce the reference code of the
paper's Simulations 1 and 2 (`tests/testthat`; scripts and results in `dev/`). Monte Carlo coverage of
the 95% sandwich intervals is 93-97% across the settings checked.
