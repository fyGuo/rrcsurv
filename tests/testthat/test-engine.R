for (model in c("cox", "err")) {
  test_that(paste(model, "log likelihood, score and information match brute force"), {
    dat <- sim_cp(n = 200, model = model, seed = 7)
    rs <- rrcsurv:::rs_build(dat$start, dat$stop, dat$event)
    pm <- pair_mats(dat, rs)
    th <- c(0.4, 0.3, -0.2)                       # (beta_x, gamma_z1, gamma_z2)
    f <- function(t) ref_ll(t, model, dat)
    if (model == "cox") {
      W <- cbind(pm$x, pm$Z); X <- matrix(0, nrow(W), 0); thw <- th
    } else {
      W <- pm$Z; X <- pm$x; thw <- c(th[2:3], th[1])
    }
    p <- rrcsurv:::pl_parts(thw, W, X, rs)
    expect_equal(p$ll, f(th), tolerance = 1e-10)
    o <- if (model == "cox") 1:3 else c(3, 1, 2)
    expect_equal(unname(p$U[o]), num_grad(f, th), tolerance = 1e-6)
    H <- t(vapply(1:3, function(k) num_grad(function(t) {
      tw <- if (model == "cox") t else c(t[2:3], t[1])
      rrcsurv:::pl_parts(tw, W, X, rs)$U[o][k]
    }, th), numeric(3)))
    expect_equal(unname(p$I[o, o]), -H, tolerance = 1e-6)
  })
}

test_that("Cox fit matches survival::coxph with Breslow ties", {
  dat <- sim_cp(n = 500, seed = 11)
  rs <- rrcsurv:::rs_build(dat$start, dat$stop, dat$event)
  pm <- pair_mats(dat, rs)
  f <- rrcsurv:::fit_risk_model("cox", pm$x, pm$Z, rs)
  ref <- survival::coxph(survival::Surv(start, stop, event) ~ x + z1 + z2, data = dat, ties = "breslow")
  expect_true(f$converged)
  expect_equal(unname(f$coef), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(unname(f$vcov), unname(vcov(ref)), tolerance = 1e-6)
  expect_equal(unname(f$loglik["final"]), ref$loglik[2], tolerance = 1e-10)
})

test_that("ERR fit solves the score equation and recovers beta", {
  dat <- sim_cp(n = 3000, model = "err", beta = 1, seed = 5)
  rs <- rrcsurv:::rs_build(dat$start, dat$stop, dat$event)
  pm <- pair_mats(dat, rs)
  f <- rrcsurv:::fit_risk_model("err", pm$x, pm$Z, rs)
  expect_true(f$converged)
  expect_false(f$boundary)
  expect_lt(max(abs(f$score)), 1e-6)
  expect_named(f$coef, c("x", "z1", "z2"))
  # true values (1, 0.5, -0.3) within 3 SE
  expect_true(all(abs(f$coef - c(1, 0.5, -0.3)) < 3 * sqrt(diag(f$vcov))))
})

test_that("ERR without exposure equals Cox on the confounders", {
  dat <- sim_cp(n = 300, seed = 2)
  rs <- rrcsurv:::rs_build(dat$start, dat$stop, dat$event)
  pm <- pair_mats(dat, rs)
  a <- rrcsurv:::pl_fit(pm$Z, matrix(0, nrow(pm$Z), 0), rs)
  b <- rrcsurv:::fit_risk_model("err", matrix(0, nrow(pm$Z), 1), pm$Z, rs)
  expect_equal(unname(a$coef), unname(b$coef[2:3]), tolerance = 1e-8)
})
