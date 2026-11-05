library(survival)

base <- function() {
  d <- sim_me("additive", "cox", n = 200, n_val = 300, seed = 31)
  d$dat$s <- TRUE
  d
}
fit <- function(d, formula = Surv(start, stop, event) ~ rrcx(C) + z1, cal = c ~ C, ...) {
  rrc(formula, data = d$dat, id = "id", validation = d$val, calibration = cal, ...)
}

test_that("arguments are checked", {
  d <- base()
  expect_error(fit(d, min_pairs = -1), "min_pairs")
  expect_error(fit(d, min_pairs = NA), "min_pairs")
  expect_error(fit(d, joint = NA), "joint must be TRUE or FALSE")
  expect_error(fit(d, naive = "yes"), "naive must be TRUE or FALSE")
  expect_error(fit(d, variance = "bootstrap", nboot = 1.5), "nboot")
  expect_error(fit(d, control = list(maxiter = 5)), "unknown control element\\(s\\): maxiter")
  expect_error(fit(d, control = list(eps = -1)), "control\\$eps")
  expect_error(fit(d, start = c(0, 0, 0)), "start must hold 2 numbers: the exposure, then z1")
  expect_error(fit(d, val_time = "nope"), "val_time must name a column")
  expect_error(rrc(Surv(start, stop, event) ~ rrcx(C), data = d$dat, id = "subject",
                   validation = d$val, calibration = c ~ C), "id must name a column of data")
  v <- d$val; v$id <- NULL
  expect_error(rrc(Surv(start, stop, event) ~ rrcx(C), data = d$dat, id = "id", validation = v,
                   calibration = c ~ C), "validation needs the id column")
  expect_error(rrc(~ rrcx(C), data = d$dat, id = "id", validation = d$val, calibration = c ~ C),
               "two-sided formula")
})

test_that("unsupported formula terms are refused", {
  d <- base()
  expect_error(fit(d, Surv(start, stop, event) ~ rrcx(C) + strata(z1)), "strata\\(\\) terms")
  expect_error(fit(d, Surv(start, stop, event) ~ rrcx(C) + cluster(id)), "cluster\\(\\) terms")
  expect_error(fit(d, Surv(start, stop, event) ~ rrcx(C) + offset(z2)), "offset\\(\\) terms")
  expect_error(fit(d, Surv(start, stop, event) ~ rrcx(C) * z1), "main effect")
  expect_error(fit(d, Surv(stop, event) ~ rrcx(C)), "Surv\\(start, stop, event\\)")
})

test_that("problems in the main study data are reported", {
  d <- base()
  d1 <- d; d1$dat$start[3] <- NA
  expect_error(fit(d1), "1 row\\(s\\) of data have a missing start, stop or event \\(rows 3\\)")
  d1 <- d; d1$dat$start[3] <- d1$dat$stop[3]
  expect_error(fit(d1), "in the Surv\\(\\) response")
  d1 <- d; d1$dat$z1[c(2, 9)] <- NA
  expect_error(fit(d1), "missing values in confounder\\(s\\) z1")
  d1 <- d; d1$dat$z3 <- 2 * d1$dat$z1
  expect_error(fit(d1, Surv(start, stop, event) ~ rrcx(C) + z1 + z3), "collinear with the others: z3")
  d1 <- d; d1$dat$k <- 1
  expect_error(fit(d1, Surv(start, stop, event) ~ rrcx(C) + k), "constant or collinear")
  d1 <- d; d1$dat$id[5] <- NA
  expect_error(fit(d1), "1 row\\(s\\) of data have a missing id")
  d1 <- d; d1$dat$C[5] <- Inf
  expect_error(fit(d1), "missing or non-finite rrcx\\(\\) exposure")
  d1 <- d; d1$dat$event <- 0
  expect_error(fit(d1), "no events")

  # overlapping rows, two events, rows after the event
  i <- d$dat$id[which(d$dat$event == 1)[1]]
  r <- which(d$dat$id == i)
  d1 <- d; d1$dat <- rbind(d1$dat, d1$dat[r[1], ])
  expect_error(fit(d1), sprintf("rows of 1 subject\\(s\\) overlap in time \\(id %s\\)", i))
  j <- d$dat$id[which(d$dat$event == 0 & duplicated(d$dat$id))[1]]
  d1 <- d; rj <- which(d1$dat$id == j); d1$dat$event[rj[1]] <- 1
  expect_error(fit(d1), "rows after their event")
  extra <- d$dat[r[length(r)], ]; extra$start <- extra$stop; extra$stop <- extra$stop + 1
  d1 <- d; d1$dat <- rbind(d1$dat, extra)
  expect_error(fit(d1), "more than one event")

  expect_warning(fit(d, Surv(start, stop, event) ~ rrcx(C, surrogate = !s) + z1),
                 "no rows are flagged as surrogate")
})

test_that("problems in the validation data are reported", {
  d <- base()
  v <- d$val; v$id[1:3] <- c(-1, -2, -3)
  expect_error(rrc(Surv(start, stop, event) ~ rrcx(C), data = d$dat, id = "id", validation = v,
                   calibration = c ~ C), "3 validation pair\\(s\\) have ids not in data \\(-1, -2, -3\\)")
  expect_error(fit(d, cal = c ~ C + z9), "calibration variables not in validation: z9")
  d1 <- d; d1$dat$grp <- ifelse(d1$dat$z1 == 1, "a", "b"); d1$val$grp <- "a"
  d1$val$grp[1:2] <- "c"
  expect_error(fit(d1, cal = c ~ C + grp), "factor levels in data must also appear in validation")
  d1 <- d; d1$val <- d1$val[1:2, ]
  expect_error(fit(d1, cal = c ~ C + z1 + z2), "too few usable validation pairs \\(2\\)")
  d1 <- d; d1$val <- d1$val[1:15, ]
  expect_warning(fit(d1), "only 15 usable validation pairs")
})

test_that("non-convergence is warned about", {
  d <- base()
  expect_warning(fit(d, control = list(maxit = 1)), "did not converge")
})

test_that("surrogate_floor raises low surrogate doses before a log calibration only", {
  d <- sim_me("multiplicative", "cox", n = 200, n_val = 300, seed = 32)
  d$dat$s <- d$dat$id %% 3 > 0                     # some rows keep their dose unchanged
  low <- which(d$dat$s)[1:20]
  d$dat$C[low] <- -0.01
  fm <- Surv(start, stop, event) ~ rrcx(C, surrogate = s) + z1
  run <- function(dat, ...) rrc(fm, data = dat, id = "id", validation = d$val,
                                calibration = log(c) ~ log(C), ...)
  expect_error(run(d$dat), "see surrogate_floor")
  expect_error(run(d$dat, surrogate_floor = 0), "surrogate_floor must be positive")
  expect_error(rrc(fm, data = d$dat, id = "id", validation = d$val, calibration = c ~ C,
                   surrogate_floor = 1e-3), "multiplicative calibration")
  f <- run(d$dat, surrogate_floor = 1e-3)
  n_low <- sum(d$dat$s & d$dat$C < 1e-3)           # the 20 negative ones plus small positive ones
  expect_gte(n_low, 20)
  expect_equal(f$n[["surrogates_floored"]], n_low)
  expect_output(print(summary(f)), sprintf("%d surrogate dose\\(s\\) in data raised to surrogate_floor = 0.001", n_low))
  # same as raising them by hand; the naive fit keeps the original doses
  d2 <- d$dat; d2$C[d2$s & d2$C < 1e-3] <- 1e-3
  g <- run(d2)
  expect_equal(coef(f), coef(g))
  expect_equal(vcov(f), vcov(g))
  expect_false(isTRUE(all.equal(f$naive$coef, g$naive$coef)))
  # non-surrogate rows are not floored (their dose is the true value, used as is)
  d3 <- d$dat; d3$s[low[1]] <- FALSE
  expect_equal(run(d3, surrogate_floor = 1e-3)$n[["surrogates_floored"]], n_low - 1)
})
