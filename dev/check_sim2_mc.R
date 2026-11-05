# Phase 5 Monte Carlo: bias, SE and coverage of rrc() at the full Simulation 2 sizes (n1 = 6000,
# n2 = 2000 with up to 8 pairs each, 60 years, scenario beta1: true beta = 1), for the correctly
# specified RRC fit of two cells at rho(c, C) = 0.6:
#   cell 2:  Cox data, additive error        -> Cox with c ~ C
#   cell 11: ERR data, multiplicative error  -> ERR with log(c) ~ log(C)
# Data come from sim2_repro's code and seeds, so replicate r here is replicate r of Simulation 2.
# run from rrcsurv/:  Rscript dev/check_sim2_mc.R [reps] [cores]
args <- commandArgs(TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 300
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1, parallel::detectCores() - 1)
Sys.setenv(SIM2_SCENARIO = "beta1")
source("dev/sim2_data.R")

one_mc <- function(ci, rep) {
  cell <- CELLS[CELLS$cell == ci, ]
  s <- sim2_draw(cell, rep)
  m <- cell$model
  t0 <- proc.time()[3]
  if (cell$error == "add") {
    d <- to_rrc(s$ex, s$sv, s$pairs); cal <- c ~ dose
  } else {
    d <- to_rrc(s$ex, s$sv, s$pairs[s$pos, ]); cal <- log(c) ~ log(dose)
  }
  f <- suppressWarnings(rrc(FM2, data = d$dat, id = "id", validation = d$val, calibration = cal,
                            model = m, min_pairs = MIN_RS_PAIRS, start = BETA[[m]]))
  se <- function(type) sqrt(vcov(f, type)[1, 1])
  data.frame(cell = ci, model = m, error = cell$error, rep = rep, truth = unname(BETA[[m]]),
             events = f$n[["events"]], b = unname(coef(f)), naive = unname(f$naive$coef[1]),
             se_sandwich = se("separate"), se_joint = se("joint"), se_main = se("main_only"),
             se_model = se("model"), converged = f$fit$converged, secs = proc.time()[3] - t0)
}

jobs <- expand.grid(rep = seq_len(B), cell = c(2, 11))
t0 <- proc.time()[3]
res <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
  tryCatch(one_mc(jobs$cell[i], jobs$rep[i]), error = function(e) {
    data.frame(cell = jobs$cell[i], rep = jobs$rep[i], error_msg = conditionMessage(e))
  })
}, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(res, function(r) "error_msg" %in% names(r), logical(1))
if (any(failed)) { cat("failed replicates:\n"); print(do.call(rbind, res[failed])) }
res <- do.call(rbind, res[!failed])
saveRDS(res, "dev/check_sim2_mc.rds")

z <- qnorm(0.975)
cat(sprintf("%d replicates per cell, %.1f min on %d cores\n", B, (proc.time()[3] - t0) / 60, cores))
for (ci in unique(res$cell)) {
  r <- res[res$cell == ci, ]
  ok <- r[r$converged, ]
  cat(sprintf("\n== cell %d: %s data, %s error, RRC %s: %d of %d converged, mean %.0f events, %.1f s per fit\n",
              ci, r$model[1], r$error[1], r$model[1], nrow(ok), nrow(r), mean(r$events), median(r$secs)))
  cat(sprintf("  truth %.3f   mean %.4f (rel bias %+.1f%%)   median %.4f   naive mean %.4f (validated years at the true dose)   MC SD %.4f\n",
              ok$truth[1], mean(ok$b), 100 * (mean(ok$b) / ok$truth[1] - 1), median(ok$b),
              mean(ok$naive), sd(ok$b)))
  tab <- t(vapply(c("se_sandwich", "se_joint", "se_main", "se_model"), function(v) {
    c(mean_se = mean(ok[[v]]), se_over_sd = mean(ok[[v]]) / sd(ok$b),
      coverage = mean(abs(ok$b - ok$truth) <= z * ok[[v]]))
  }, numeric(3)))
  print(round(tab, 3))
}
