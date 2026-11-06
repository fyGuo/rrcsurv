# Cross-check: the rrcsurv engine on Simulation 1 data, in person-time format, reproduces
# sim1_repro's pl_fit() (the code behind Table 1) for the naive and true-exposure fits.
# run from rrcsurv/:  Rscript dev/check_engine_vs_sim1.R
pkgload::load_all(".", quiet = TRUE)
sim <- new.env(); sim$PROJ <- normalizePath("../measurementERROR_COX_ERR_RCC/sim1_repro")
for (f in c("params.R", "exposure.R", "survival.R", "legacy/fit.R")) sys.source(file.path("../measurementERROR_COX_ERR_RCC/sim1_repro/R", f), envir = sim)

# sim1 matrices -> person-time rows: subject i, interval (j-1, j] carries cumulative exposure x_{j-1}
to_person_time <- function(time, D, x) {
  n <- length(time); last <- ceiling(time)
  id <- rep(seq_len(n), last); j <- sequence(last)
  data.frame(id = id, start = j - 1, stop = pmin(j, time[id]),
             event = as.integer(D[id] == 1 & j == last[id]), x = x[cbind(id, j)])
}

with(sim, {
  set.seed(1)
  L <- cum_lag_matrix(T_STAR)
  for (cell in list(list(model = "cox", error = "add"), list(model = "err", error = "mult"))) {
    cl <- CELLS[CELLS$model == cell$model & CELLS$error == cell$error & CELLS$cor == "ind", ]
    ex <- gen_exposure(N_MAIN, T_STAR, cell$error, COR$ind, expo_params(cell$error, cl$sd, cl$rho))
    x <- ex$c %*% L; X <- ex$C %*% L
    sv <- gen_survival(x, cell$model, BETA[[cell$model]], cl$nu, THETA)
    ev <- make_events(sv$time, sv$D, sv$k)
    for (nm in c("true", "naive")) {
      M <- if (nm == "true") x else X
      old <- pl_fit(z_long(M, ev), ev$eid, z_case(M, ev), cell$model, start = 0)
      pt <- to_person_time(sv$time, sv$D, M)
      rs <- rs_build(pt$start, pt$stop, pt$event)
      new <- fit_risk_model(cell$model, matrix(pt$x[rs$row], dimnames = list(NULL, "x")), NULL, rs)
      cat(sprintf("%s/%s %-5s  sim1 beta %.10f se %.10f | rrcsurv beta %.10f se %.10f | diff %.1e\n",
                  cell$model, cell$error, nm, old$beta, old$se, new$coef, sqrt(new$vcov),
                  abs(old$beta - new$coef)))
    }
  }
})
