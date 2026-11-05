library(survival)
cells2 <- readRDS(test_path("fixtures", "sim2_cells.rds"))

# Simulation 2 design: 60 years of follow-up, up to 8 validation pairs per subject; reference values
# from sim2_repro's rrc_fit() (dev/make_sim2_fixture.R)
for (cell in cells2) {
  test_that(paste("Simulation 2 cell", cell$cell, "(", cell$model, cell$error, "): rrc() reproduces rrc_fit()"), {
    fA <- rrc(Surv(start, stop, event) ~ rrcx(dose_A, surrogate = surr_A), data = cell$main, id = "id",
              validation = transform(cell$validation, dose_A = dose), calibration = c ~ dose_A,
              model = cell$model, min_pairs = cell$min_pairs, start = cell$start, naive = FALSE)
    # log scale: non-positive pairs are dropped and, under additive data, non-positive surrogates floored
    vM <- transform(cell$validation, dose_M = dose)
    run_M <- function() {
      rrc(Surv(start, stop, event) ~ rrcx(dose_M, surrogate = surr_M), data = cell$main, id = "id",
          validation = vM, calibration = log(c) ~ log(dose_M), model = cell$model,
          min_pairs = cell$min_pairs, start = cell$start, naive = FALSE,
          surrogate_floor = if (cell$error == "add") 1e-8)
    }
    if (any(!cell$validation$pos)) {
      expect_warning(fM <- run_M(), sprintf("%d of %d validation pairs", sum(!cell$validation$pos),
                                            nrow(cell$validation)))
      expect_gt(fM$n[["surrogates_floored"]], 0)
    } else {
      fM <- run_M()
    }
    for (x in list(list(f = fA, ref = cell$ref$A), list(f = fM, ref = cell$ref$M))) {
      expect_true(x$f$fit$converged)
      expect_equal(nrow(x$f$groups), x$ref$n_groups)
      expect_equal(unname(coef(x$f)), x$ref$beta, tolerance = 1e-6)
      expect_equal(sqrt(drop(vcov(x$f))), x$ref$se, tolerance = 1e-6)
      expect_equal(sqrt(drop(vcov(x$f, "main_only"))), x$ref$se_beta_only, tolerance = 1e-6)
      expect_equal(sqrt(drop(vcov(x$f, "model"))), x$ref$se_model, tolerance = 1e-6)
    }
  })
}
