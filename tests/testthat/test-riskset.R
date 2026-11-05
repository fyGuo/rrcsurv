test_that("risk sets match a brute-force construction", {
  dat <- sim_cp(n = 150, seed = 3)
  rs <- rcrsurv_rs <- rrcsurv:::rs_build(dat$start, dat$stop, dat$event)
  expect_equal(rs$time, sort(unique(dat$stop[dat$event == 1])))
  for (k in seq_along(rs$time)) {
    t <- rs$time[k]
    expect_setequal(rs$row[rs$uid == k], which(dat$start < t & dat$stop >= t))
    expect_equal(rs$d[k], sum(dat$event == 1 & dat$stop == t))
  }
  # each case points at its own row at its own event time
  cases <- which(dat$event == 1)
  expect_setequal(rs$row[rs$case], cases)
  expect_equal(rs$time[rs$uid[rs$case]], dat$stop[rs$row[rs$case]])
})

test_that("rs_build rejects bad input", {
  expect_error(rrcsurv:::rs_build(c(0, 1), c(1, 1), c(0, 1)), "start < stop")
  expect_error(rrcsurv:::rs_build(c(0, 1), c(1, 2), c(0, 0)), "no events")
})
