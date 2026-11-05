library(survival)
cells <- readRDS(test_path("fixtures", "sim1_cells.rds"))

for (cell in cells) {
  test_that(paste(cell$model, cell$error, "sandwich reproduces the Simulation 1 rrc_fit() SEs"), {
    cf <- if (cell$error == "add") c ~ dose else log(c) ~ log(dose)
    f <- rrc(Surv(start, stop, event) ~ rrcx(dose, surrogate = surr), data = cell$main, id = "id",
             validation = cell$validation, calibration = cf, model = cell$model,
             min_pairs = cell$min_pairs)
    expect_equal(sqrt(drop(vcov(f))), cell$ref$se, tolerance = 1e-7)
    expect_equal(sqrt(drop(f$variances$main_only)), cell$ref$se_beta_only, tolerance = 1e-7)
    expect_equal(sqrt(drop(f$variances$model)), cell$ref$se_model, tolerance = 1e-7)
    m <- cell$ref_merge
    f2 <- rrc(Surv(start, stop, event) ~ rrcx(dose, surrogate = surr), data = cell$main, id = "id",
              validation = cell$validation, calibration = cf, model = cell$model,
              min_pairs = m$min_pairs)
    expect_equal(sqrt(drop(vcov(f2))), m$se, tolerance = 1e-7)
  })
}

combos <- expand.grid(model = c("cox", "err"), scale = c("additive", "multiplicative"),
                      stringsAsFactors = FALSE)
for (i in seq_len(nrow(combos))) {
  model <- combos$model[i]; scale <- combos$scale[i]
  test_that(paste(model, scale, "with confounders: U_i sum to the score and U* = dU/dalpha"), {
    d <- sim_me(scale, model, seed = 10 + i)
    cf <- if (scale == "additive") c ~ C + z1 + z2 else log(c) ~ log(C) + z1
    f <- rrc(Surv(start, stop, event) ~ rrcx(C) + z1 + z2, data = d$dat, id = "id",
             validation = d$val, calibration = cf, model = model, min_pairs = 150)
    expect_true(f$fit$converged)
    it <- f$internal
    sp <- rrcsurv:::rrc_score_parts(f$fit, it$xrrc, it$rs, it$sid, it$grp, it$fits, it$k, it$first,
                                    it$surr, it$dose, it$Mm, f$scale)
    expect_equal(unname(colSums(sp$Ui)), unname(f$fit$score), tolerance = 1e-8)
    expect_gt(length(it$fits), 1)
    for (g in unique(c(1, length(it$fits)))) {
      fg <- it$fits[[g]]
      par <- c(fg$alpha, if (scale == "multiplicative") fg$s2)
      num <- vapply(seq_along(par), function(j) {
        h <- 1e-5 * max(1, abs(par[j])); e <- replace(numeric(length(par)), j, h)
        (score_at_alpha(f, g, par + e) - score_at_alpha(f, g, par - e)) / (2 * h)
      }, numeric(length(f$coefficients)))
      expect_equal(unname(sp$Ustar[[g]]), unname(num), tolerance = 1e-6)
    }
    V <- vcov(f)
    expect_equal(dimnames(V), list(names(coef(f)), names(coef(f))))
    expect_true(all(diag(V) >= diag(f$variances$main_only) - 1e-12))
  })
}

test_that("joint = TRUE changes the sandwich only through the internal-validation cross terms", {
  d <- sim_me("additive", "cox", seed = 3)
  a <- rrc(Surv(start, stop, event) ~ rrcx(C) + z1, data = d$dat, id = "id", validation = d$val,
           calibration = c ~ C + z1, min_pairs = 150)
  b <- rrc(Surv(start, stop, event) ~ rrcx(C) + z1, data = d$dat, id = "id", validation = d$val,
           calibration = c ~ C + z1, min_pairs = 150, joint = TRUE)
  expect_equal(coef(a), coef(b))
  expect_false(isTRUE(all.equal(vcov(a), vcov(b))))
  expect_equal(sqrt(diag(vcov(a))), sqrt(diag(vcov(b))), tolerance = 0.1)
  n <- rrc(Surv(start, stop, event) ~ rrcx(C) + z1, data = d$dat, id = "id", validation = d$val,
           calibration = c ~ C + z1, min_pairs = 150, variance = "none")
  expect_true(all(is.na(vcov(n))))
})

test_that("bootstrap SE agrees with the sandwich, internal and external validation", {
  skip_on_cran()
  d <- sim_me("multiplicative", "cox", n = 600, n_val = 700, seed = 21)
  fm <- Surv(start, stop, event) ~ rrcx(C) + z1 + z2
  set.seed(1)
  f <- rrc(fm, data = d$dat, id = "id", validation = d$val, calibration = log(c) ~ log(C) + z1,
           min_pairs = 150, variance = "bootstrap", nboot = 80)
  expect_equal(f$boot$n_failed, 0)
  expect_identical(vcov(f), f$variances$bootstrap)
  ratio <- sqrt(diag(f$variances$bootstrap) / diag(f$variances$sandwich))
  expect_true(all(ratio > 0.7 & ratio < 1.4))

  # external: validation subjects are not in the main study and carry their own follow-up time
  ev <- d$val
  ev$T <- tapply(d$dat$stop, d$dat$id, max)[as.character(ev$id)]
  ev$id <- ev$id + 1e5
  set.seed(2)
  e <- rrc(fm, data = d$dat, id = "id", validation = ev, calibration = log(c) ~ log(C) + z1,
           min_pairs = 150, val_time = "T", variance = "bootstrap", nboot = 80)
  expect_equal(coef(e), coef(f))
  expect_equal(e$variances$sandwich, f$variances$sandwich)   # same pairs, same risk sets
  ratio <- sqrt(diag(e$variances$bootstrap) / diag(e$variances$sandwich))
  expect_true(all(ratio > 0.7 & ratio < 1.4))
})
