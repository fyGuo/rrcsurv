# Shared by dev/check_sim2.R and dev/check_sim2_mc.R: load sim2_repro (the code behind Table 2) and
# rrcsurv, regenerate one Simulation 2 replicate and convert it to rrc() input.
# Source from rrcsurv/ after setting SIM2_SCENARIO.

# load.R uses source(), which evaluates in the global environment; rrc() uses its own namespace
PROJ <- normalizePath("../sim2_repro")
local({ owd <- setwd(PROJ); on.exit(setwd(owd)); source("R/load.R") })
pkgload::load_all(".", quiet = TRUE)
library(survival)

# The data of one_rep(cell, rep): one_rep() draws exposure, survival, then validation after set.seed
sim2_draw <- function(cell, rep, base_seed = 20261102) {
  set.seed(base_seed + 100000 * cell$cell + rep)
  L <- cum_lag_matrix(T_STAR)
  ex <- gen_exposure2(N_MAIN, T_STAR, cell$error, expo_params(cell$error, cell$rho))
  x <- ex$c %*% L
  sv <- gen_survival(x, cell$model, BETA[[cell$model]], cell$nu, THETA)
  pairs <- sample_validation2(sv$time, N_VAL_SUBJ, MAX_PER_ID, T_STAR)
  # the log calibration uses only pairs with c, C > 0 (one_rep())
  pos <- ex$c[cbind(pairs$id, pairs$l)] > 0 & ex$C[cbind(pairs$id, pairs$l)] > 0
  list(ex = ex, sv = sv, pairs = pairs, pos = pos)
}

# Person-time data for rrc(): subject i, year j = interval (j - 1, j], ending at the event or
# censoring time. Validated (subject, year) pairs carry the true dose and surrogate = FALSE, as in
# rrc_fit(), which uses the observed c in validated years. The cumulative exposure for (k, k + 1]
# is the dose of years 1..k, which is rrcx(lag = 0).
to_rrc <- function(ex, sv, pairs) {
  n <- nrow(ex$C)
  last <- ceiling(sv$time)
  id <- rep(seq_len(n), last); j <- sequence(last)
  val <- paste(id, j) %in% paste(pairs$id, pairs$l)
  dat <- data.frame(id = id, start = j - 1, stop = pmin(j, sv$time[id]),
                    event = as.integer(sv$D[id] == 1 & j == last[id]),
                    dose = ifelse(val, ex$c[cbind(id, j)], ex$C[cbind(id, j)]), surr = !val)
  v <- data.frame(id = pairs$id, c = ex$c[cbind(pairs$id, pairs$l)], dose = ex$C[cbind(pairs$id, pairs$l)])
  list(dat = dat, val = v)
}

FM2 <- Surv(start, stop, event) ~ rrcx(dose, surrogate = surr)
