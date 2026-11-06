library(survival)

# Validation risk sets with dose times: at event time t, the pairs of subjects at risk at t
# (T >= t) whose dose entered the cumulative exposure by t (u + lag < t).

fit_dt <- function(d, scale = "additive", model = "cox", lag = 0, ...) {
  cf <- if (scale == "additive") c ~ C + z1 else log(c) ~ log(C) + z1
  fm <- if (lag > 0) Surv(start, stop, event) ~ rrcx(C, lag = lag) + z1 else
    Surv(start, stop, event) ~ rrcx(C) + z1
  rrc(fm, data = d$dat, id = "id", validation = d$val, calibration = cf, model = model, ...)
}

# brute force: the pairs of each risk-set group, from the definition
brute_sets <- function(f, Tv, u, lag) {
  times <- f$internal$rs$time; grp <- f$internal$grp
  lapply(seq_len(max(grp)), function(g) {
    tg <- times[grp == g]
    which(vapply(seq_along(Tv), function(j) any(u[j] + lag < tg & Tv[j] >= tg), logical(1)))
  })
}

test_that("dose times long before every event time reproduce the default risk sets", {
  d <- sim_me("additive", "cox", seed = 21)
  a <- fit_dt(d, min_pairs = 150)
  d$val$u_old <- -1e6
  b <- fit_dt(d, min_pairs = 150, dose_time = "u_old")
  expect_equal(coef(b), coef(a))
  expect_equal(vcov(b), vcov(a))
  expect_equal(b$groups, a$groups)
  expect_equal(unname(b$n["pairs_unused"]), unname(a$n["pairs_unused"]))
})

for (lag in c(0, 2)) {
  test_that(paste("risk sets with dose times match the definition, lag", lag), {
    d <- sim_me("additive", "cox", seed = 22)
    f <- fit_dt(d, lag = lag, min_pairs = 100, dose_time = "u")
    expect_true(f$fit$converged)
    Tsub <- tapply(d$dat$stop, d$dat$id, max)
    Tv <- unname(Tsub[as.character(d$val$id)])
    sets <- brute_sets(f, Tv, d$val$u, lag)
    expect_gt(length(sets), 1)
    for (g in seq_along(sets)) {
      S <- f$internal$fits[[g]]$S
      expect_equal(S, sets[[g]])
      expect_gt(length(S), 100)
      ref <- lm(c ~ C + z1, data = d$val[S, ])
      expect_equal(unname(f$internal$fits[[g]]$alpha), unname(coef(ref)), tolerance = 1e-10)
    }
    used <- sort(unique(unlist(sets)))
    expect_equal(unname(f$n["pairs_unused"]), sum(!vapply(seq_along(Tv), function(j)
      any(d$val$u[j] + lag < f$internal$rs$time & Tv[j] >= f$internal$rs$time), logical(1))))
    expect_true(length(used) <= nrow(d$val) - f$n["pairs_unused"])
    # no dose is received before the first event time, so its risk set is empty and it shares
    # a group with later event times
    expect_length(brute_sets(list(internal = list(rs = list(time = f$internal$rs$time[1]), grp = 1L)),
                             Tv, d$val$u, lag)[[1]], 0)
    expect_gt(sum(f$internal$grp == 1), 1)
    a <- fit_dt(d, lag = lag, min_pairs = 100)
    expect_false(isTRUE(all.equal(coef(f), coef(a))))
  })
}

combos <- expand.grid(model = c("cox", "err"), scale = c("additive", "multiplicative"),
                      stringsAsFactors = FALSE)
for (i in seq_len(nrow(combos))) {
  model <- combos$model[i]; scale <- combos$scale[i]
  test_that(paste(model, scale, "with dose times: U* = dU/dalpha in every group"), {
    d <- sim_me(scale, model, seed = 30 + i)
    f <- fit_dt(d, scale, model, min_pairs = 120, dose_time = "u")
    expect_true(f$fit$converged)
    it <- f$internal
    sp <- rrcsurv:::rrc_score_parts(f$fit, it$xrrc, it$rs, it$sid, it$grp, it$fits, it$k, it$first,
                                    it$surr, it$dose, it$Mm, f$scale)
    expect_equal(unname(colSums(sp$Ui)), unname(f$fit$score), tolerance = 1e-8)
    for (g in seq_along(it$fits)) {
      fg <- it$fits[[g]]
      par <- c(fg$alpha, if (scale == "multiplicative") fg$s2)
      num <- vapply(seq_along(par), function(j) {
        h <- 1e-5 * max(1, abs(par[j])); e <- replace(numeric(length(par)), j, h)
        (score_at_alpha(f, g, par + e) - score_at_alpha(f, g, par - e)) / (2 * h)
      }, numeric(length(f$coefficients)))
      expect_equal(unname(sp$Ustar[[g]]), unname(num), tolerance = 1e-6)
    }
    expect_true(all(diag(vcov(f)) >= diag(f$variances$main_only) - 1e-12))
  })
}

test_that("external validation study: dose times combine with val_time", {
  d <- sim_me("additive", "cox", seed = 23)
  Tsub <- tapply(d$dat$stop, d$dat$id, max)
  d$val$T <- unname(Tsub[as.character(d$val$id)])
  d$val$id <- -d$val$id
  f <- fit_dt(d, min_pairs = 100, val_time = "T", dose_time = "u")
  expect_false(f$internal_validation)
  sets <- brute_sets(f, d$val$T, d$val$u, 0)
  for (g in seq_along(sets)) expect_equal(f$internal$fits[[g]]$S, sets[[g]])
})

test_that("the bootstrap keeps the dose-time rule", {
  d <- sim_me("additive", "cox", seed = 24, n = 300, n_val = 400)
  d$val$u_old <- -1e6
  set.seed(5); a <- fit_dt(d, min_pairs = 80, variance = "bootstrap", nboot = 15)
  set.seed(5); b <- fit_dt(d, min_pairs = 80, variance = "bootstrap", nboot = 15, dose_time = "u_old")
  set.seed(5); u <- fit_dt(d, min_pairs = 80, variance = "bootstrap", nboot = 15, dose_time = "u")
  expect_equal(b$boot$coef, a$boot$coef)
  expect_false(isTRUE(all.equal(u$boot$coef, a$boot$coef)))
})

test_that("summary reports the risk-set rule", {
  d <- sim_me("additive", "cox", seed = 25)
  f <- fit_dt(d, min_pairs = 100, dose_time = "u")
  out <- capture.output(summary(f))
  expect_true(any(grepl("with dose received before t (dose_time = \"u\")", out, fixed = TRUE)))
  expect_true(any(grepl("pair(s) in no risk set", out, fixed = TRUE)))
  out0 <- capture.output(summary(fit_dt(d, min_pairs = 100)))
  expect_true(any(grepl("(any dose time)", out0, fixed = TRUE)))
})

test_that("dose_time input problems give clear messages", {
  d <- sim_me("additive", "cox", seed = 26)
  expect_error(fit_dt(d, dose_time = "nope"), "dose_time must name a column")
  expect_error(fit_dt(d, dose_time = c("u", "u")), "dose_time must name a column")
  d2 <- d; d2$val$u <- as.character(d2$val$u)
  expect_error(fit_dt(d2, dose_time = "u"), "must be numeric")
  d3 <- d; d3$val$u[1:3] <- NA
  expect_error(fit_dt(d3, dose_time = "u"), "3 validation pair\\(s\\) have a missing dose_time")
  d4 <- d; d4$val$u <- 1e6                        # every dose after the last event
  expect_error(fit_dt(d4, dose_time = "u"), "too few usable validation pairs \\(0; 600 more are in no risk set\\)")
})
