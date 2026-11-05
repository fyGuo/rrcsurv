# Build tests/testthat/fixtures/sim1_cells.rds: small Simulation 1 data sets (one per risk model x
# error model, independent correlation) in person-time format, with the estimates of sim1_repro's
# rrc_fit() (the code behind Table 1) as reference values.
# run from rrcsurv/:  Rscript dev/make_sim1_fixture.R
sim <- new.env(); sim$PROJ <- normalizePath("../sim1_repro")
for (f in c("params.R", "exposure.R", "survival.R", "fit.R", "validation.R", "rrc.R")) {
  sys.source(file.path(sim$PROJ, "R", f), envir = sim)
}
n1 <- 600; n2 <- 200; n_pair <- 750

out <- with(sim, lapply(1:4, function(i) {
  cell <- CELLS[CELLS$cor == "ind", ][i, ]
  set.seed(1000 + i)
  L <- cum_lag_matrix(T_STAR)
  ex <- gen_exposure(n1, T_STAR, cell$error, COR$ind, expo_params(cell$error, cell$sd, cell$rho))
  x <- ex$c %*% L
  sv <- gen_survival(x, cell$model, BETA[[cell$model]], cell$nu, THETA)
  ev <- make_events(sv$time, sv$D, sv$k)
  pairs <- sample_validation(sv$time, n2, n_pair, MAX_PER_ID, T_STAR)
  ref <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, cell$error, cell$model, L, MIN_RS_PAIRS, start = 0)
  ref2 <- rrc_fit(pairs, ex$C, ex$c, sv$time, ev, cell$error, cell$model, L, 690, start = 0)  # forces merging

  # person-time rows: subject i, year j = interval (j - 1, j]; validated years carry the true dose
  last <- ceiling(sv$time)
  id <- rep(seq_len(n1), last); j <- sequence(last)
  val <- paste(id, j) %in% paste(pairs$id, pairs$l)
  main <- data.frame(id = id, start = j - 1, stop = pmin(j, sv$time[id]),
                     event = as.integer(sv$D[id] == 1 & j == last[id]),
                     dose = ifelse(val, ex$c[cbind(id, j)], ex$C[cbind(id, j)]), surr = !val)
  validation <- data.frame(id = pairs$id, c = ex$c[cbind(pairs$id, pairs$l)],
                           dose = ex$C[cbind(pairs$id, pairs$l)])
  list(model = cell$model, error = cell$error, main = main, validation = validation,
       min_pairs = MIN_RS_PAIRS,
       ref = list(beta = ref$beta, se = ref$se, se_beta_only = ref$se_beta_only, se_model = ref$se_model,
                  n_groups = ref$n_groups, converged = ref$converged),
       ref_merge = list(min_pairs = 690, beta = ref2$beta, se = ref2$se, n_groups = ref2$n_groups))
}))
saveRDS(out, "tests/testthat/fixtures/sim1_cells.rds", compress = "xz")
for (o in out) cat(sprintf("%s/%s: %d rows, %d events, beta %.6f, se %.6f, groups %d; merged: beta %.6f groups %d\n", o$model, o$error,
                           nrow(o$main), sum(o$main$event), o$ref$beta, o$ref$se, o$ref$n_groups, o$ref_merge$beta, o$ref_merge$n_groups))
