library(survival)

fit_me <- function(model, scale = "additive", ...) {
  d <- sim_me(scale, model, seed = 21)
  cf <- if (scale == "additive") c ~ C + z1 else log(c) ~ log(C) + z1
  rrc(Surv(start, stop, event) ~ rrcx(C) + z1 + z2, data = d$dat, id = "id",
      validation = d$val, calibration = cf, model = model, min_pairs = 150, ...)
}

test_that("vcov() gives each stored variance by type", {
  f <- fit_me("cox")
  expect_identical(vcov(f), f$vcov)
  expect_identical(vcov(f, "sandwich"), f$variances$sandwich)
  expect_identical(vcov(f, "separate"), f$variances$sandwich_separate)
  expect_identical(vcov(f, "joint"), f$variances$sandwich_joint)
  expect_identical(vcov(f, "model"), f$variances$model)
  expect_identical(vcov(f, "main"), f$variances$main_only)        # partial matching
  expect_error(vcov(f, "bootstrap"), "variance = \"bootstrap\"")
  expect_error(vcov(f, "nonsense"))
  fj <- fit_me("cox", joint = TRUE)
  expect_identical(vcov(fj), fj$variances$sandwich_joint)
  expect_output(print(fj), "joint sandwich")
})

test_that("confint() is the Wald interval for any level, parm and type", {
  for (model in c("cox", "err")) {
    f <- fit_me(model, "multiplicative")
    se <- sqrt(diag(vcov(f)))
    ci <- confint(f)
    expect_equal(unname(ci[, 1]), unname(coef(f) - qnorm(0.975) * se))
    expect_equal(unname(ci[, 2]), unname(coef(f) + qnorm(0.975) * se))
    expect_equal(colnames(ci), c("2.5 %", "97.5 %"))
    ci2 <- confint(f, "C", level = 0.9, type = "model")
    expect_equal(dim(ci2), c(1L, 2L))
    expect_equal(unname(ci2[1, 2] - ci2[1, 1]), 2 * qnorm(0.95) * sqrt(f$variances$model[1, 1]))
    expect_equal(confint(f, 2:3), ci[2:3, ])
    expect_error(confint(f, "w"), "parm")
    expect_error(confint(f, level = 1), "below 1")
  }
})

test_that("summary() tables: Wald tests, HR for Cox, ERR per unit for the ERR exposure", {
  f <- fit_me("cox")
  s <- summary(f)
  tab <- s$coefficients
  expect_equal(tab[, "coef"], coef(f))
  expect_equal(tab[, "z"], coef(f) / sqrt(diag(vcov(f))))
  expect_equal(tab[, "Pr(>|z|)"], 2 * pnorm(-abs(tab[, "z"])))
  expect_equal(unname(s$conf.int), unname(exp(cbind(coef(f), confint(f)))))
  expect_equal(s$naive[, "coef"], f$naive$coef)
  expect_output(print(s), "exp\\(coef\\)")
  expect_output(print(s), "Calibration parameters across risk-set groups")

  e <- fit_me("err")
  se <- summary(e, conf.level = 0.9)
  ci <- confint(e, level = 0.9)
  expect_equal(unname(se$conf.int["C", ]), unname(c(coef(e)["C"], ci["C", ])))
  expect_equal(unname(se$conf.int[c("z1", "z2"), ]),
               unname(exp(cbind(coef(e), ci)[c("z1", "z2"), ])))
  expect_equal(colnames(se$conf.int), c("estimate", "lower .90", "upper .90"))
  expect_output(print(se), "excess relative rate per unit")
  # summary on another variance
  expect_equal(summary(e, type = "model")$coefficients[, "se(coef)"], sqrt(diag(vcov(e, "model"))))
  expect_equal(nobs(e), unname(e$n["events"]))
})

test_that("an external validation study has no joint sandwich; variance = 'none' has no SEs", {
  d <- sim_me("additive", "cox", seed = 22)
  Tsub <- tapply(d$dat$stop, d$dat$id, max)
  v <- d$val
  v$T <- as.numeric(Tsub[as.character(v$id)])
  v$id <- -seq_len(nrow(v))
  f <- rrc(Surv(start, stop, event) ~ rrcx(C) + z1, data = d$dat, id = "id", validation = v,
           calibration = c ~ C + z1, val_time = "T", min_pairs = 150)
  expect_false(f$internal_validation)
  expect_error(vcov(f, "joint"), "internal validation")
  expect_output(print(summary(f)), "External validation study")
  expect_warning(rrc(Surv(start, stop, event) ~ rrcx(C) + z1, data = d$dat, id = "id", validation = v,
                     calibration = c ~ C + z1, val_time = "T", min_pairs = 150, joint = TRUE),
                 "internal validation study only")

  n <- rrc(Surv(start, stop, event) ~ rrcx(C) + z1, data = d$dat, id = "id", validation = d$val,
           calibration = c ~ C + z1, variance = "none", naive = FALSE)
  expect_true(all(is.na(vcov(n))))
  expect_error(vcov(n, "sandwich"), "variance = \"none\"")
  expect_output(print(n), "not computed")
  expect_output(print(summary(n)), "Calibration parameters")
})
