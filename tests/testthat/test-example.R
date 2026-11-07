test_that("Surv() is re-exported and found without survival on the search path", {
  expect_identical(rrcsurv::Surv, survival::Surv)
  # formula environment that sees neither survival nor rrcsurv (as on a parallel worker)
  fm <- stats::as.formula("Surv(start, stop, event) ~ rrcx(dose_true) + smoke",
                          env = new.env(parent = baseenv()))
  f <- cumfit(fm, data = rrc_main, id = "id")
  ref <- cumfit(survival::Surv(start, stop, event) ~ rrcx(dose_true) + smoke, data = rrc_main, id = "id")
  expect_equal(coef(f), coef(ref))
})

test_that("the example data match their documentation", {
  expect_equal(nrow(rrc_main), 22947)
  expect_equal(sum(rrc_main$event), 482)
  expect_equal(length(unique(rrc_main$id)), 1500)
  expect_equal(nrow(rrc_validation), 1159)
  expect_equal(length(unique(rrc_validation$id)), 400)
  # validated years carry the true dose; validation's year is the stop of the pair's row
  m <- match(paste(rrc_validation$id, rrc_validation$year), paste(rrc_main$id, rrc_main$stop))
  expect_false(anyNA(m))
  expect_true(all(!rrc_main$surrogate[m]))
  expect_equal(sum(!rrc_main$surrogate), nrow(rrc_validation))
  expect_equal(rrc_main$dose[m], rrc_validation$c)
  expect_equal(rrc_main$dose_true[m], rrc_validation$c)
})

test_that("RRC on the example data recovers the true-dose fit", {
  fit <- rrc(Surv(start, stop, event) ~ rrcx(dose, surrogate = surrogate) + smoke, data = rrc_main,
             id = "id", validation = rrc_validation, calibration = log(c) ~ log(dose) + smoke,
             dose_time = "year")
  truth <- cumfit(Surv(start, stop, event) ~ rrcx(dose_true) + smoke, data = rrc_main, id = "id")
  expect_true(fit$fit$converged)
  expect_lt(abs(coef(fit)[["dose"]] - coef(truth)[["dose_true"]]), 2 * sqrt(vcov(fit)[1, 1]))
  expect_lt(fit$naive$coef[["dose"]], 0.6 * coef(fit)[["dose"]])
})
