library(survival)

test_that("cumfit() Cox equals coxph: estimate, model-based and robust SE, with a lag", {
  dat <- sim_cp(n = 300, seed = 41)
  dat <- dat[order(dat$id, dat$start), ]
  for (lag in c(0, 2)) {
    dat$xl <- vapply(seq_len(nrow(dat)), function(j) {
      sum(dat$dose[dat$id == dat$id[j] & dat$stop <= dat$start[j] - lag])
    }, 0)
    f <- cumfit(Surv(start, stop, event) ~ rrcx(dose, lag = lag) + z1 + z2, data = dat, id = "id")
    ref <- coxph(Surv(start, stop, event) ~ xl + z1 + z2, data = dat, ties = "breslow")
    rob <- coxph(Surv(start, stop, event) ~ xl + z1 + z2, data = dat, ties = "breslow", cluster = id)
    expect_true(f$fit$converged)
    expect_equal(unname(coef(f)), unname(coef(ref)), tolerance = 1e-8)
    expect_equal(unname(vcov(f)), unname(vcov(ref)), tolerance = 1e-6)
    expect_equal(unname(vcov(f, "robust")), unname(vcov(rob)), tolerance = 1e-6)
    expect_named(coef(f), c("dose", "z1", "z2"))
  }
})

test_that("cumfit() ERR maximises the partial likelihood and equals rrc()'s naive fit", {
  dat <- sim_cp(n = 400, model = "err", seed = 42)
  f <- cumfit(Surv(start, stop, event) ~ rrcx(dose) + z1 + z2, data = dat, id = "id", model = "err")
  expect_true(f$fit$converged)
  dat$x <- ave(dat$dose, dat$id, FUN = cumsum) - dat$dose
  th <- unname(coef(f))
  expect_lt(max(abs(num_grad(function(t) ref_ll(t, "err", dat), th))), 1e-5)
  H <- vapply(seq_along(th), function(k) {
    e <- replace(numeric(3), k, 1e-4)
    (num_grad(function(t) ref_ll(t, "err", dat), th + e) - num_grad(function(t) ref_ll(t, "err", dat), th - e)) / 2e-4
  }, numeric(3))
  expect_equal(unname(vcov(f)), solve(-H), tolerance = 1e-4)

  v <- data.frame(id = dat$id[1:200], c = dat$dose[1:200], dose = dat$dose[1:200])
  r <- rrc(Surv(start, stop, event) ~ rrcx(dose) + z1 + z2, data = dat, id = "id", validation = v,
           calibration = c ~ dose, model = "err", variance = "none")
  expect_equal(coef(f), r$naive$coef, tolerance = 1e-10)
})

test_that("cumfit() methods", {
  dat <- sim_cp(n = 300, model = "err", seed = 43)
  f <- cumfit(Surv(start, stop, event) ~ rrcx(dose) + z1, data = dat, id = "id", model = "err")
  ci <- confint(f, level = 0.9, type = "robust")
  expect_equal(unname(ci[, 2] - ci[, 1]), unname(2 * qnorm(0.95) * sqrt(diag(vcov(f, "robust")))))
  s <- summary(f)
  expect_equal(s$coefficients[, "se(coef)"], sqrt(diag(vcov(f))))
  expect_equal(unname(s$conf.int["z1", ]), unname(exp(c(coef(f)["z1"], confint(f)["z1", ]))))
  expect_equal(unname(s$conf.int["dose", 1]), unname(coef(f)["dose"]))
  expect_output(print(s), "no measurement error correction")
  expect_output(print(f), "Linear ERR model")
  expect_equal(nobs(f), sum(dat$event))
  expect_error(vcov(f, "sandwich"))
  expect_error(cumfit(Surv(start, stop, event) ~ dose, data = dat, id = "id"), "exactly one rrcx")
  expect_warning(cumfit(Surv(start, stop, event) ~ rrcx(dose), data = dat, id = "id",
                        control = list(maxit = 1)), "did not converge")
})
