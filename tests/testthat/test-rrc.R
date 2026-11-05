library(survival)
cells <- readRDS(test_path("fixtures", "sim1_cells.rds"))

for (cell in cells) {
  test_that(paste(cell$model, cell$error, "RRC reproduces the Simulation 1 rrc_fit() estimate"), {
    cf <- if (cell$error == "add") c ~ dose else log(c) ~ log(dose)
    for (ref in list(c(cell$ref, min_pairs = cell$min_pairs), cell$ref_merge)) {
      f <- rrc(Surv(start, stop, event) ~ rrcx(dose, surrogate = surr), data = cell$main, id = "id",
               validation = cell$validation, calibration = cf, model = cell$model,
               min_pairs = ref$min_pairs)
      expect_true(f$fit$converged)
      expect_equal(unname(coef(f)), ref$beta, tolerance = 1e-7)
      expect_equal(nrow(f$groups), ref$n_groups)
    }
  })
}

# person-time data with a validation subset whose true dose equals the surrogate
cp_with_validation <- function(seed = 4, n = 400) {
  dat <- sim_cp(n = n, seed = seed)
  dat$smoke <- factor(ifelse(dat$z1 == 1, "yes", "no"))
  set.seed(seed)
  v <- dat[sample(nrow(dat), 500), c("id", "dose", "z1", "z2", "smoke")]
  list(dat = dat, val = v)
}

test_that("with c = C, RRC equals the naive fit, and the naive fit equals coxph", {
  d <- cp_with_validation()
  d$val$c <- d$val$dose
  f <- rrc(Surv(start, stop, event) ~ rrcx(dose) + smoke + z2, data = d$dat, id = "id",
           validation = d$val, calibration = c ~ dose + z2 + smoke, model = "cox", min_pairs = 50)
  expect_equal(coef(f), f$naive$coef, tolerance = 1e-8)
  expect_true(all(abs(as.matrix(f$groups[, c("dose", "z2", "smokeyes")]) - rep(c(1, 0, 0), each = nrow(f$groups))) < 1e-10))
  d$dat$xprev <- ave(d$dat$dose, d$dat$id, FUN = cumsum) - d$dat$dose
  ref <- coxph(Surv(start, stop, event) ~ xprev + smoke + z2, data = d$dat, ties = "breslow")
  expect_equal(unname(f$naive$coef), unname(coef(ref)), tolerance = 1e-7)
  expect_named(coef(f), c("dose", "smokeyes", "z2"))
})

test_that("ERR with c = C under the multiplicative model also equals naive", {
  d <- cp_with_validation(seed = 9)
  d$val$c <- d$val$dose
  f <- suppressWarnings(rrc(Surv(start, stop, event) ~ rrcx(dose) + z1, data = d$dat, id = "id",
                            validation = d$val, calibration = log(c) ~ log(dose), model = "err"))
  # log c = log C exactly: alpha = (0, 1), sigma2 = 0, so E(c | C) = C
  expect_equal(coef(f), f$naive$coef, tolerance = 1e-7)
})

test_that("cumulative exposure honours the lag", {
  dat <- sim_cp(n = 60, seed = 2)
  o <- order(dat$id, dat$start); dat <- dat[o, ]
  sid <- match(dat$id, unique(dat$id))
  for (lag in c(0, 1, 2.5)) {
    k <- rrcsurv:::prefix_index(sid, dat$start, dat$stop, lag)
    x <- rrcsurv:::cum_from_prefix(dat$dose, k, match(sid, sid))
    ref <- vapply(seq_len(nrow(dat)), function(j) {
      sum(dat$dose[dat$id == dat$id[j] & dat$stop <= dat$start[j] - lag])
    }, 0)
    expect_equal(x, ref)
  }
})

test_that("risk-set groups merge backwards until each has more than min_pairs pairs", {
  expect_equal(rrcsurv:::rrc_groups(c(30, 25, 15, 10, 5), 20), c(1L, 2L, 2L, 2L, 2L))
  expect_equal(rrcsurv:::rrc_groups(c(30, 25, 21), 20), 1:3)
  expect_equal(rrcsurv:::rrc_groups(c(30, 25, 21), Inf), c(1L, 1L, 1L))
})

test_that("input problems give clear messages", {
  d <- cp_with_validation()
  d$val$c <- d$val$dose
  d$val$c[1:5] <- -1
  expect_warning(rrc(Surv(start, stop, event) ~ rrcx(dose), data = d$dat, id = "id", validation = d$val,
                     calibration = log(c) ~ log(dose)), "5 of 500 validation pairs")
  v2 <- d$val; names(v2)[names(v2) == "dose"] <- "C"
  expect_error(rrc(Surv(start, stop, event) ~ rrcx(dose), data = d$dat, id = "id", validation = v2,
                   calibration = c ~ C), "same name as the rrcx")
  expect_error(rrc(Surv(start, stop, event) ~ dose, data = d$dat, id = "id", validation = d$val,
                   calibration = c ~ dose), "exactly one rrcx")
  expect_error(rrc(Surv(start, stop, event) ~ rrcx(dose), data = d$dat, id = "id", validation = d$val,
                   calibration = sqrt(c) ~ dose), "left side")
  v3 <- d$val; v3$id <- v3$id + 10000
  expect_error(rrc(Surv(start, stop, event) ~ rrcx(dose), data = d$dat, id = "id", validation = v3,
                   calibration = c ~ dose), "external validation")
})
