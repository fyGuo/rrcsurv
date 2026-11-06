# Build tests/testthat/fixtures/sim2_cells.rds: two small Simulation 2 data sets (sim2_repro, the code
# behind Table 2) in person-time format, with sim2_repro's rrc_fit() estimates and SEs as reference:
#   cell 2  (Cox, additive error, rho 0.6): cox_A, and cox_M on the log scale with the floored surrogate
#   cell 11 (ERR, multiplicative error, rho 0.6): err_A and err_M
# Several validation pairs per subject, 60 years of follow-up. Sizes are smaller than in the
# simulation (n1 = 1500, n2 = 400) to keep the fixture and the tests small.
# run from rrcsurv/:  Rscript dev/make_sim2_fixture.R
Sys.setenv(SIM2_SCENARIO = "beta1")
PROJ <- normalizePath("../measurementERROR_COX_ERR_RCC/sim2_repro")
local({ owd <- setwd(PROJ); on.exit(setwd(owd)); source("R/load.R") })
n1 <- 1500; n2 <- 400

out <- lapply(c(2, 11), function(ci) {
  cell <- CELLS[CELLS$cell == ci, ]
  set.seed(5000 + ci)
  L <- cum_lag_matrix(T_STAR)
  ex <- gen_exposure2(n1, T_STAR, cell$error, expo_params(cell$error, cell$rho))
  x <- ex$c %*% L
  sv <- gen_survival(x, cell$model, BETA[[cell$model]], cell$nu, THETA)
  ev <- make_events(sv$time, sv$D, sv$k)
  pairs <- sample_validation2(sv$time, n2, MAX_PER_ID, T_STAR)
  pos <- ex$c[cbind(pairs$id, pairs$l)] > 0 & ex$C[cbind(pairs$id, pairs$l)] > 0
  m <- cell$model
  refA <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, "add", m, L, MIN_RS_PAIRS, start = BETA[[m]])
  refM <- rrc_fit(pairs[pos, ], pmax(ex$C, 1e-8), ex$c, sv$time, ev, "mult", m, L, MIN_RS_PAIRS,
                  start = BETA[[m]])

  # person-time rows; validated years carry the true dose. For the log calibration only pairs
  # with c, C > 0 count as validated (surr_M), as in sim2_repro's one_rep().
  last <- ceiling(sv$time)
  id <- rep(seq_len(n1), last); j <- sequence(last)
  key <- paste(id, j)
  vA <- key %in% paste(pairs$id, pairs$l)
  vM <- key %in% paste(pairs$id[pos], pairs$l[pos])
  main <- data.frame(id = id, start = j - 1, stop = pmin(j, sv$time[id]),
                     event = as.integer(sv$D[id] == 1 & j == last[id]),
                     C = ex$C[cbind(id, j)], c = ex$c[cbind(id, j)], surr_A = !vA, surr_M = !vM)
  main$dose_A <- ifelse(vA, main$c, main$C)
  main$dose_M <- ifelse(vM, main$c, main$C)
  main$C <- main$c <- NULL
  validation <- data.frame(id = pairs$id, c = ex$c[cbind(pairs$id, pairs$l)],
                           dose = ex$C[cbind(pairs$id, pairs$l)], pos = pos)
  list(cell = ci, model = m, error = cell$error, main = main, validation = validation,
       min_pairs = MIN_RS_PAIRS, start = unname(BETA[[m]]),
       ref = list(A = refA[c("beta", "se", "se_beta_only", "se_model", "n_groups", "converged")],
                  M = refM[c("beta", "se", "se_beta_only", "se_model", "n_groups", "converged")]))
})
str(lapply(out, function(o) c(rows = nrow(o$main), events = sum(o$main$event),
                              pairs = nrow(o$validation), unlist(o$ref))))
saveRDS(out, "tests/testthat/fixtures/sim2_cells.rds", compress = "xz")
cat(file.size("tests/testthat/fixtures/sim2_cells.rds") / 1e6, "MB\n")
