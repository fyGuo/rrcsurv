# Monte Carlo check of the RRC sandwich variance with confounders (phase 3).
# Continuous event times, confounders z1 (binary, also shifts the dose) and z2, yearly doses
# until a random end of exposure (zero, error-free dose after it),
# a calibration model that is exactly right, internal validation of some subject-years.
# For each model x error cell: mean SE (sandwich as in the paper, joint sandwich, main-study part
# only, model-based) vs the Monte Carlo SD of the estimates, and 95% CI coverage.
# run from rrcsurv/:  Rscript dev/check_variance_mc.R [R]
pkgload::load_all(".", quiet = TRUE)
library(survival); library(parallel)
args <- as.integer(commandArgs(TRUE))
R_REPS <- if (length(args) >= 1) args[1] else 500

LAMBDA <- list(cox = 0.018, err = 0.022)             # about 15% with an event
TRUE_PAR <- list(cox = c(dose = log(1.5), z1 = 0.4, z2 = -0.3),
                 err = c(dose = 0.5, z1 = 0.4, z2 = -0.3))

gen <- function(model, scale, n = 3000, T = 20, lambda = LAMBDA[[model]], n_val_subj = 500, per_subj = 3) {
  b <- TRUE_PAR[[model]]
  z1 <- rbinom(n, 1, 0.5); z2 <- rnorm(n)
  id <- rep(seq_len(n), each = T); j <- rep(seq_len(T), n); zz1 <- z1[id]
  expo <- j <= sample(0:T, n, replace = TRUE)[id]                    # exposed years: until a random year
  if (scale == "additive") {
    C <- rgamma(n * T, shape = 4, rate = 20 / (1 + 0.5 * zz1))       # mean 0.2, 0.3 if z1
    c <- 0.02 + 0.6 * C + 0.04 * zz1 + rnorm(n * T, 0, 0.05)
  } else {
    C <- exp(rnorm(n * T, log(0.2) + 0.2 * zz1, 0.6))
    c <- exp(-0.5 + 0.8 * log(C) + 0.15 * zz1 + rnorm(n * T, 0, 0.4))
  }
  C[!expo] <- 0; c[!expo] <- 0                                       # unexposed years: known zero dose
  # true cumulative dose before year j, and the hazard in year j
  cm <- matrix(c, T); x <- rbind(0, apply(cm, 2, cumsum)[-T, , drop = FALSE])
  lp <- b["z1"] * z1 + b["z2"] * z2
  rr <- if (model == "cox") exp(b["dose"] * x + rep(lp, each = T)) else
    pmax(1 + b["dose"] * x, 1e-8) * exp(rep(lp, each = T))
  h0 <- lambda * (seq_len(T) / T)^3                                  # rising baseline, constant within a year
  H <- apply(h0 * rr, 2, cumsum)                                     # cumulative hazard at year ends
  E <- rexp(n)
  cens <- pmin(T, rexp(n, 0.03))
  ev_year <- apply(H >= rep(E, each = T), 2, function(v) if (any(v)) which(v)[1] else NA)
  Hprev <- ifelse(is.na(ev_year), NA, H[cbind(pmax(ev_year - 1, 1), seq_len(n))] * (ev_year > 1))
  t_ev <- ifelse(is.na(ev_year), Inf,
                 ev_year - 1 + (E - Hprev) / (h0[ev_year] * rr[cbind(ev_year, seq_len(n))]))
  time <- pmin(t_ev, cens); D <- as.integer(t_ev <= cens)
  keep <- j - 1 < time[id]
  main <- data.frame(id = id, start = j - 1, stop = pmin(j, time[id]),
                     event = as.integer(D[id] == 1 & j >= time[id]), C = C, c = c,
                     z1 = zz1, z2 = z2[id], expo = expo)[keep, ]
  # internal validation: up to per_subj observed years of n_val_subj random subjects
  vs <- sample(n, n_val_subj)
  vr <- unlist(lapply(vs, function(i) {
    r <- which(main$id == i & main$expo); r[sample.int(length(r), min(per_subj, length(r)))]
  }))
  main$surr <- main$expo; main$surr[vr] <- FALSE
  val <- main[vr, c("id", "c", "C", "z1", "z2")]
  main$dose <- ifelse(main$surr, main$C, main$c)                    # validated years: true dose
  list(main = main[, c("id", "start", "stop", "event", "dose", "surr", "z1", "z2")],
       val = transform(val, dose = C))
}

one <- function(rep, model, scale) {
  set.seed(20261105 + 1000 * rep + 100 * (model == "err") + 10 * (scale == "multiplicative"))
  d <- gen(model, scale)
  cf <- if (scale == "additive") c ~ dose + z1 else log(c) ~ log(dose) + z1
  f <- tryCatch(rrc(Surv(start, stop, event) ~ rrcx(dose, surrogate = surr) + z1 + z2, data = d$main,
                    id = "id", validation = d$val, calibration = cf, model = model, naive = FALSE),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  v <- f$variances
  data.frame(rep = rep, par = names(coef(f)), est = coef(f), se = sqrt(diag(v$sandwich)),
             se_joint = sqrt(diag(v$sandwich_joint)), se_main = sqrt(diag(v$main_only)),
             se_model = sqrt(diag(v$model)), conv = f$fit$converged, events = f$n["events"],
             groups = nrow(f$groups))
}

cells <- expand.grid(scale = c("additive", "multiplicative"), model = c("cox", "err"), stringsAsFactors = FALSE)
all <- list()
for (i in seq_len(nrow(cells))) {
  m <- cells$model[i]; s <- cells$scale[i]
  t0 <- Sys.time()
  res <- do.call(rbind, mclapply(seq_len(R_REPS), function(r) one(r, m, s),
                                 mc.cores = max(1, detectCores() - 1)))
  n_fail <- R_REPS - length(unique(res$rep[res$conv]))
  res <- res[res$conv, ]
  res$model <- m; res$scale <- s; res$truth <- TRUE_PAR[[m]][res$par]
  all[[i]] <- res
  cat(sprintf("\n== %s / %s: %d converged reps (%d failed or not converged), mean %.0f events, %.1f min\n",
              m, s, length(unique(res$rep)), n_fail, mean(res$events),
              as.numeric(Sys.time() - t0, units = "mins")))
  smry <- do.call(rbind, lapply(split(res, res$par), function(r) {
    cover <- function(se) mean(abs(r$est - r$truth) <= 1.96 * se)
    data.frame(par = r$par[1], truth = r$truth[1], rel_bias = mean(r$est) / r$truth[1] - 1,
               mc_sd = sd(r$est), se_sandwich = mean(r$se), se_joint = mean(r$se_joint),
               se_main = mean(r$se_main), se_model = mean(r$se_model),
               cover_sandwich = cover(r$se), cover_joint = cover(r$se_joint),
               cover_main = cover(r$se_main), cover_model = cover(r$se_model))
  }))
  print(smry, digits = 3, row.names = FALSE, width = 200)
}
saveRDS(do.call(rbind, all), "dev/check_variance_mc.rds")
