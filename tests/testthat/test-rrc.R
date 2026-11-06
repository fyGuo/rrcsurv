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

# nested risk sets (no dose times) holding np[e] pairs at event time e: every pair starts at 1
nested <- function(np) {
  hi <- rep(seq_along(np), np - c(np[-1], 0))
  list(lo = rep(1L, length(hi)), hi = hi)
}
groups_of <- function(r, E, min_pairs) rrcsurv:::rrc_groups(r$lo, r$hi, E, min_pairs)

test_that("risk-set groups merge backwards until each has more than min_pairs pairs", {
  expect_equal(groups_of(nested(c(30, 25, 15, 10, 5)), 5, 20), c(1L, 2L, 2L, 2L, 2L))
  expect_equal(groups_of(nested(c(30, 25, 21)), 3, 20), 1:3)
  expect_equal(groups_of(nested(c(30, 25, 21)), 3, Inf), c(1L, 1L, 1L))
})

test_that("non-nested risk sets: a group holds the union, and a small earliest group joins the next", {
  # 5 pairs only at event time 1, 30 only at 2, 25 only at 3
  r <- list(lo = rep(1:3, c(5, 30, 25)), hi = rep(1:3, c(5, 30, 25)))
  expect_equal(groups_of(r, 3, 20), c(1L, 1L, 2L))
  # 15 pairs at 3, 10 spanning 2..3, 30 at 1: event time 3 alone has 25 > 20; 2 has 10 and
  # merges with 1 (10 + 30 = 40), never counting a pair twice
  r <- list(lo = c(rep(3, 15), rep(2, 10), rep(1, 30)), hi = c(rep(3, 15), rep(3, 10), rep(1, 30)))
  expect_equal(groups_of(r, 3, 20), c(1L, 1L, 2L))
  expect_equal(groups_of(r, 3, 24), c(1L, 1L, 2L))
  expect_equal(groups_of(r, 3, 25), c(1L, 1L, 1L))
  # pairs in no risk set (lo > hi) are ignored
  r2 <- list(lo = c(r$lo, rep(3, 50)), hi = c(r$hi, rep(2, 50)))
  expect_equal(groups_of(r2, 3, 20), groups_of(r, 3, 20))
})

test_that("rrc_groups matches a brute-force backward merge on random pair ranges", {
  brute <- function(lo, hi, E, m) {
    grp <- integer(E); g <- 0L; b <- NA
    for (e in E:1) {
      if (is.na(b)) { g <- g + 1L; b <- e }
      grp[e] <- g
      if (sum(lo <= b & hi >= e & lo <= hi) > m) b <- NA
    }
    if (!is.na(b) && g > 1L) { grp[grp == g] <- g - 1L; g <- g - 1L }
    g + 1L - grp
  }
  set.seed(7)
  for (rep in 1:200) {
    E <- sample(1:12, 1); P <- sample(0:60, 1)
    lo <- sample(1:(E + 1), P, replace = TRUE); hi <- pmin(E, lo + sample(-1:E, P, replace = TRUE))
    m <- sample(c(0, 3, 8, 20), 1)
    expect_equal(rrcsurv:::rrc_groups(lo, hi, E, m), brute(lo, hi, E, m))
  }
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
